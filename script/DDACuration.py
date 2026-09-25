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

"""DDACuration.py

Automatic curation of MassBank record files produced from DDA (Data-Dependent
Acquisition) MS2 data.

This script assigns a confidence level (Level 3-SP / Level 3-CL / Level 4 /
Level 5, following the MSI/MassBank annotation confidence scheme) to each
averaged MS2 spectrum by inspecting the candidate compounds, molecular
formulas, and compound classes proposed by upstream in-silico annotation
tools (SIRIUS, MS-FINDER, MetFrag, CFM-ID, GNPS spectral similarity, etc.).
Based on the assigned level, the corresponding MassBank record is updated
with curated compound information (name, formula, links), a curated record
title, fragment-ion formula annotations, and a human-readable confidence
comment, and is then written out as a finalized MassBank record file.

Typical usage (from the command line)::

    python DDACuration.py HSA001 HSA002

If no dataset is given, every dataset found under the data folder is
processed.
"""

import ddatoolbox
import glob, copy, pickle, csv, itertools, sys, datetime

from statistics import stdev
from pathlib import Path
import os, re
import argparse

def updateCHfield(mb2, pdcompound):
    """Updates the chemical information (CH$) fields in a MassBank record.

    Args:
        mb2 (MassBankRecord): The MassBank record object to be updated.
        pdcompound (dict): A dictionary containing PubChem compound data,
            including 'IUPACName', 'MolecularFormula', 'IsomericSMILES',
            'InChI', 'ExactMass', 'CID', and 'InChIKey'.

    Returns:
        MassBankRecord: The updated MassBank record object.
    """
    if not pdcompound == {}:
        mb2.clear_anything("CH$NAME:", "*")
        IUPACName = "No IUPAC Name"
        if 'IUPACName' in pdcompound:
            IUPACName =pdcompound['IUPACName']
        mb2.add("CH$NAME:", IUPACName)
        mb2.add("CH$COMPOUND_CLASS:", "N/A")
        mb2.add("CH$FORMULA:", pdcompound['MolecularFormula'])
        mb2.add("CH$SMILES:", pdcompound['IsomericSMILES'])
        mb2.add("CH$IUPAC:", pdcompound['InChI'])
        mb2.add("CH$EXACT_MASS:", pdcompound['ExactMass'])
        mb2.add2("CH$LINK:", "PUBCHEM", "CID:"+str(pdcompound['CID']))
        mb2.add2("CH$LINK:", "INCHIKEY", pdcompound['InChIKey'])
    return(mb2)



def fragmentformulasearch(precmf, mb2, logfilename):
    """Searches for molecular formulas for fragment ions based on the precursor formula.

    Args:
        precmf (str): Molecular formula of the precursor ion.
        mb2 (MassBankRecord): The MassBank record object containing spectral data.
        logfilename (str): Path to the log file for recording search results.

    Returns:
        tuple: A tuple containing:
            - temp_hash (dict): Updated annotation data for the fragments.
            - spectra_comp_fragment (dict): A dictionary of fragment m/z and intensities
              that matched a formula.
    """
    mzthreshold_ms2 = float(mb2.get_anything("COMMENT:", "MS2_M/Z_AVERAGED_THRES"))
    precmz, precint = mb2.get_precmz()
    polarity = mb2.get_polarity()


    ionform = mb2.get_anything("MS$FOCUSED_ION:", "ION_TYPE")

    result = ddatoolbox.mfsearch.elementnumber(precmf)
    Climit = result["C"]
    Hlimit = result["H"]
    Nlimit = result["N"]
    Olimit = result["O"]
    Plimit = result["P"]
    Slimit = result["S"]


    temp_hash = mb2.get_annotation_from_items()

    temp_fragmentmz = sorted(temp_hash.keys())
    spectra_comp_fragment = {}

    for fragmentmz in temp_fragmentmz:

        temp_hash[fragmentmz]["mzPAF"] = temp_hash[fragmentmz]["annotation"]

        threshold_temp = mzthreshold_ms2/1000000 * fragmentmz

        fragmentintensity = int(temp_hash[fragmentmz]["intensity"])
        fragmentmzstdev = 0 #float(temp_hash[fragmentmz]["mzsd"])
        fragmentnumber = float(temp_hash[fragmentmz]["number"])
        correlation = float(temp_hash[fragmentmz]["correlation"])

        if float(fragmentmz) > precmz + threshold_temp * 2:
            logtemp =["  ", "{:6f} {:3d} {:8f} ".format(fragmentmz,fragmentintensity, correlation), "above prec"]
            #print(logtemp)
            ddatoolbox.ddatoolbox.writeLog(logfilename, logtemp)

            continue
        #
        # Ignore isotope (13C) candidates
        #
        if ddatoolbox.mfsearch.check13C(fragmentmz, temp_hash, threshold = threshold_temp * 2):
            logtemp =["  ", "{:6f} {:3d} {:8f} ".format(fragmentmz,fragmentintensity,  correlation),  "  ignored"]
            #print(logtemp)
            ddatoolbox.ddatoolbox.writeLog(logfilename, logtemp)
            continue
        #
        # sodium and pottasium adducts
        #
        if "+Na" in ionform:
            fragmentmz_temp = fragmentmz - 22.98977 # sodium atom
        elif "+K" in ionform:
            fragmentmz_temp = fragmentmz - 38.96371 # pottasium atom
        else:
            fragmentmz_temp = fragmentmz

        #
        # Fragment formula search
        #

        fragment_result, fragment_result_list = ddatoolbox.mfsearch.mfsearch(fragmentmz_temp, threshold = float(mzthreshold_ms2)*2, ihdtype = "both", \
        Climit = Climit , Hlimit =Hlimit, Nlimit = Nlimit, Olimit = Olimit, Plimit = Plimit, Slimit = Slimit, polarity = polarity)

        if len(fragment_result_list) == 0 and fragmentmz_temp != fragmentmz:
            fragmentmz_temp = fragmentmz
            fragment_result, fragment_result_list = ddatoolbox.mfsearch.mfsearch(fragmentmz_temp, threshold = float(mzthreshold_ms2)*2, ihdtype = "both", \
            Climit = Climit , Hlimit =Hlimit, Nlimit = Nlimit, Olimit = Olimit, Plimit = Plimit, Slimit = Slimit, polarity = polarity)

        if len(fragment_result_list) > 0:
            #
            # compariton of precmz
            #
            #print(fragmentmz, fragment_result_list)
            for mftemp in fragment_result_list[0:1]:
                exactprecmz = ddatoolbox.mfsearch.calcmw(mftemp, polarity = polarity)
                error = (fragmentmz_temp - exactprecmz)/exactprecmz*1000000
                diff = ddatoolbox.mfsearch.diffformula(precmf, mftemp)


                if abs(error) < float(mzthreshold_ms2):
                    spectra_comp_fragment[fragmentmz] = fragmentintensity

                #
                # sodium and pottasium adducts
                #

                if diff == "":
                    difftext = "p"
                else:
                    difftext = "p-"+diff

                if fragmentmz_temp != fragmentmz:
                    if "+Na" in ionform:
                        mftemp += "Na1"
                    elif "+K" in ionform:
                        mftemp += "K1"
                else:
                    if "+Na" in ionform:
                        difftext += "Na1"
                    elif "+K" in ionform:
                        difftext += "K1"

                paf = "f{{{:s}}}/{:.2f}ppm,{:s}".format(mftemp,  error, difftext)

                temp_hash[fragmentmz]["annotation"] = paf
                temp_hash[fragmentmz]["mzPAF"] = paf
                logtemp=["  ", "{:6f} {:3d} {:8f} {:s}".format(fragmentmz,fragmentintensity, correlation, paf )]
                #print(logtemp)
                ddatoolbox.ddatoolbox.writeLog(logfilename, logtemp)

        else:
            logtemp=["--", fragmentmz, "{:3d} {:6f}".format(fragmentintensity, correlation)]
            #print(logtemp)
            ddatoolbox.ddatoolbox.writeLog(logfilename, logtemp)
    return(temp_hash, spectra_comp_fragment)




def modifyclass3(mb2, hitrecord, leveltext, pubchem_checked, InChIkeydb, logfilename, unrealistic, spec_option = ""):
    """Modifies MassBank records for Level 3 annotations using external database metadata.

    Args:
        mb2 (MassBankRecord): The MassBank record object.
        hitrecord (dict): Information about the identified hit (inchikey, hmdb, etc.).
        classtext (str): The confidence level text.
        pubchem_checked (dict): Cache of previously fetched PubChem data.
        InChIkeydb (dict): Local database of InChIKeys and compound metadata.
        logfilename (str): Path to the log file.
        spec_option (str, optional): Option to filter spectra (e.g., "hitonly"). Defaults to "".

    Returns:
        MassBankRecord: The updated MassBank record object.
    """
    #
    # Updating  CH: field by inchikey compound
    #
    accession = mb2.get_accession() #id
    #
    # Update HMDB ID
    #
    if hitrecord["hmdb"] != "":
        mb2.add2("CH$LINK:", "HMDB", hitrecord["hmdb"])
    #
    # gathering PubChem data
    #
    inchikey = hitrecord["inchikey"]
    inchikey_long = "-".join([inchikey, "UHFFFAOYSA", "N"])
    #
    # gathering namecandidates
    #
    #anno = getannodata(mb2, target = "compound")
    anno = mb2.get_proposed(target = "compound")



    namecandidates = {}
    for method, temp in anno["tg"][inchikey].items():
        namecandidates[method] = temp["name"]

    if inchikey_long in pubchem_checked:
        pdcompound = pubchem_checked[inchikey_long]
    else:
        pdcompound = ddatoolbox.resources.get_pubchemdatafrominchikey(inchikey_long)
        if pdcompound == {}:
            pdcompound = ddatoolbox.resources.get_pubchemdatafrominchikey(inchikey)
        pubchem_checked[inchikey_long] = pdcompound

    if not pdcompound == {}:
        #
        # Update of CH$NAME:
        #
        mb2 = updateCHfield(mb2, pdcompound)
        #
        # Compound datanase of MS-FINDER would be more curated than pubchem
        #
        if inchikey in InChIkeydb:
            InChIkeytemp = InChIkeydb[inchikey]
            #"Title", "Exact mass", "ChEBI","HMDB","KNApSAcK","LipidMAPS", "InChIKeyshort", "InChIkey"
            if not InChIkeytemp["Title"] == "N/A":
                namecandidates["MS-FINDER:MS-FINDER_DB"] = InChIkeytemp["Title"]
            if not InChIkeytemp["ChEBI"] == "N/A":
                templist = InChIkeytemp["ChEBI"].split(";")
                mb2.add2("CH$LINK:", "CHEBI", templist[0])
            if not InChIkeytemp["HMDB"] == "N/A":
                templist = InChIkeytemp["HMDB"].split(";")
                mb2.add2("CH$LINK:", "HMDB", templist[0])
            if not InChIkeytemp["LipidMAPS"] == "N/A":
                templist = InChIkeytemp["LipidMAPS"].split(";")
                mb2.add2("CH$LINK:", "LIPIDMAPS", templist[0])
        #
        # Updating CH$NAME by compond name in search results
        #
        # Collect up to three distinct name candidates, in priority order:
        #   name  - PubChem-derived name (from CH$NAME:, skipping synonym entries with "=")
        #   name2 - name from the MS-FINDER internal structure database (curated)
        #   name3 - shortest name proposed by MetFrag/SIRIUS confidence tools
        names = mb2.record[("CH$NAME:", "*")]["data"]
        name = ""
        for nametemp in names:
            if not nametemp == "nan":
                if not "=" in nametemp:
                    name = nametemp
                    break

        name2 = ""

        if "MS-FINDER:MS-FINDER_DB" in namecandidates:
            if not namecandidates["MS-FINDER:MS-FINDER_DB"] == "nan":
                name2 = namecandidates["MS-FINDER:MS-FINDER_DB"]

        name3 = ""
        for method in ["MetFrag:MS-FINDER_DB", "MetFrag:PubChemLite", 'SIRIUS5_Confidence', 'SIRIUS6_Confidence']:
            if method in namecandidates:
                if namecandidates[method] == "nan":
                    continue
                if name3 == "":
                    name3 = namecandidates[method]
                if len(namecandidates[method]) < len(name3):
                    # Prefer the shortest available name (usually the least decorated/systematic one)
                    name3 = namecandidates[method]
                    #print(method, name2, namecandidates[method])


        if ";" in name:
                name = name.split(";")[0]
        name = clearname(name)
        if ";" in name2:
                name2 = name2.split(";")[0]
        name2 = clearname(name2)
        if ";" in name3:
                name3 = name3.split(";")[0]
        name3 = clearname(name3)


        mb2.clear_anything("CH$NAME:", "*")
        # Rebuild CH$NAME: with the unique, non-empty candidates (name2, then name3, then name)
        namestemp = []
        if not name2 == "":
            if not name2 == "nan":
                namestemp.append(name2)
        if not name3 == name2:
            if not name3 == "":
                if not name3 == "nan":
                    namestemp.append(name3)
        if not name == name2 and not name == name3:
            if not name == "":
                if not name == "nan":
                    namestemp.append(name)
        if len(namestemp) == 0:
            namestemp.append("No name")


        mb2.record[("CH$NAME:", "*")]["data"] = namestemp

        #
        # Updating updateRecordtitile"
        #
        #level = ddatoolbox.paths.getsettings("RESOURCES", "level3sp_def")
        annotationlevel = "Confidence {}".format(leveltext)

        topname = namestemp[0]

        if inchikey in unrealistic:
           topname = "WARNING: Unrealistic structure. "+topname


        #mb2 = updateRecordtitile(mb2, namestemp[0], annotationlevel)

        mb2.add_curatedrecordtitle(topname, annotationlevel)

        #
        # Addition formula of each fragments by Peak Annotation Format
        #
        #ionformpdcompound = hitrecord["precursor_ionform"]
        ionformpdcompound = anno["data"][inchikey]["ionform"]
        mf = pdcompound['MolecularFormula']

        precmf = ddatoolbox.mfsearch.generateprecmf(mf, ionformpdcompound)

        mb2.add2("MS$FOCUSED_ION:","ION_TYPE", ionformpdcompound)

        annotation, spectra_comp_fragment = fragmentformulasearch(precmf, mb2, logfilename)

        mb2.add_annotation_from_items(annotation, ["int", "m/z_nocal", "m/z_std(ppm)", "int_std", "cor", "n", "annotation"])

        if spec_option == "hitonly":
            if len(spectra_comp_fragment) > 0:
                mb2.add_spectra(spectra_comp_fragment)



    return(mb2)



def clearname(name):
    """Cleans compound names by removing stereochemical prefixes and notations using regex.

    Args:
        name (str): The raw compound name.

    Returns:
        str: The cleaned compound name.
    """
    #
    # Strip parenthesized R/S stereo-descriptors, e.g. "(2R)-", "(2R,3S)-"
    #
    name = re.sub(r"\([0-9]*[SR]\)\-", "", name)
    name = re.sub(r"\([0-9]*[SR],[0-9]*[SR]\)\-", "", name)
    name = re.sub(r"\([0-9]*[SR],[0-9]*[SR],[0-9]*[SR]\)\-", "", name)

    #
    # Strip parenthesized E/Z and R/S locant descriptors, e.g. "(9Z,12Z)-"
    #
    name = re.sub(r"\([EZ0-9\,]+\)\-?", "", name)
    name = re.sub(r"\([RS0-9\,]+\)\-?", "", name)

    #
    # Strip bare (non-parenthesized) R/S descriptors, e.g. " R-", "R-", "R,S-"
    #
    name = re.sub(r" [RS,]+\-", "", name)
    name = re.sub(r"^[RS]\-", "", name)
    name = re.sub(r"^[RS,]+\-", "", name)
    name = re.sub(r"\-[RS,]+\-", "-", name)

    #
    # Strip bare E/Z descriptors
    #
    name = re.sub(r" [EZ,]+\-", "", name)
    name = re.sub(r"^[EZ]\-", "", name)
    name = re.sub(r"^[EZ,]+\-", "", name)
    name = re.sub(r"\-[EZ,]+\-", "-", name)

    #
    # Strip D/L descriptors (amino acid / sugar stereochemistry)
    #
    name = re.sub(r" [DL,]+\-", "", name)
    name = re.sub(r"^[DL]\-", "", name)
    name = re.sub(r"^[DL,]+\-", "", name)
    name = re.sub(r"\-[DL,]+\-", "-", name)

    #
    # Strip lowercase d/l descriptors
    #
    name = re.sub(r" [d,l]+\-", "", name)
    name = re.sub(r"^[dl]\-", "", name)
    name = re.sub(r"^[d,l]+\-", "", name)
    name = re.sub(r"\-[d,l]+\-", "-", name)

    #
    # Strip optical-rotation descriptors, e.g. "(+)-", "(-)-"
    #
    name = re.sub(r"^\([\+\-]\)\-", "", name)
    name = re.sub(r"\-?\([\+\-]\)\-", "-", name)
    return(name)


#########################################################################



def checklevel5(mb2):
    """Generates the annotation text for Level 5 confidence (Exact mass of interest).

    Args:
        mb2 (MassBankRecord): The MassBank record object.

    Returns:
        str: Formatted confidence string for Level 5.
    """
    # COMMENT: CONFIDENCE Exact mass of interest (Level 5); m/z 205.143381;  Observed in MS1 m/z 205.143381; Charge state Unknown;
    precmz, sd = mb2.get_precmz()
    # Get texts
    level = ddatoolbox.paths.getsettings("RESOURCES", "level5_def")
    level_exp = ddatoolbox.paths.getsettings("RESOURCES", "level5_exp")
    level_exp2 = ddatoolbox.paths.getsettings("RESOURCES", "level5_exp2")
    #"L5 (Exact mass of interest)"
    leveltext = "{:s} ({:s})".format(level, level_exp)


    inventory = []
    inventory.append(leveltext)
    #inventory.append("prec m/z {:.5f}".format(precmz))
    inventory.append(level_exp2.format(precmz))
    annotationtext = "; ".join(inventory)
    return(annotationtext)

def checklevel4(mb2, method = "qvalue"):
    """Checks if the record meets the criteria for Level 4 confidence (Unequivocal precursor ion formula).

    A formula candidate qualifies for Level 4 when it is ranked #1 in the
    aggregated ranking and its Bayesian q-value-like score is below 0.01.

    Args:
        mb2 (MassBankRecord): The MassBank record object.
        method (str, optional): Reserved for future use (not referenced
            inside this function; the q-value-like score is always used).
            Defaults to "qvalue".

    Returns:
        tuple: A tuple containing:
            - annotationtext (str): Formatted confidence string for Level 4
              (empty string if the record does not qualify).
            - hitrecord (dict): Dictionary with keys "formula" and "ionform"
              of the winning candidate (empty strings if no hit).
            - annotationtext2 (str): Formatted evidence-inventory string
              ("FM (Formula/mass evidence) ...") describing which tools
              ranked the candidate formula #1.
    """

    #COMMENT: CONFIDENCE Unequivocal precursor ion formula (Level 4);C10H21O4:[M+H]+; SevenGoldenRules (Ranking 1), MS-FINDER (Ranking 1),
    #SIRIUS6_Sirius (Ranking 1), SIRIUS6_Zodiac (Ranking 1);
    hitrecord = {"formula":"", "ionform":""}

    #anno = getannodata(mb2, target = "formula")
    anno = mb2.get_proposed(target = "formula")
    if not 1 in anno["rank"]:
        return("", hitrecord, "")



    #
    # check level 4 - 2 Top in aggregated ranking with q-value < 0.01
    #

    qlikevalue = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue")

    tophits = anno["rank"][1]
    hitformula = ""
    usedmethods = ['Aggregated_ranking','SIRIUS6_Sirius', 'SIRIUS5_Sirius', 'MS-FINDER']

    #if 'q-value_like' in tophits:
    if qlikevalue in tophits:

        topformula = tophits[qlikevalue]['tg']
        topscore = tophits[qlikevalue]['score']
        if topscore < 0.01:
            hitformula = topformula
        else:
            hitformula = ""
    inventory = []
    inventory2 = []
    if hitformula != "":
        diff = anno["data"][hitformula]["name"]
        ionform = anno["data"][hitformula]["ionform"]
        score = tophits[qlikevalue]['score']

        # Get texts
        level = ddatoolbox.paths.getsettings("RESOURCES", "level4_def")
        level_exp = ddatoolbox.paths.getsettings("RESOURCES", "level4_exp")
        level_exp2 = ddatoolbox.paths.getsettings("RESOURCES", "level4_exp2")
        #"L5 (Exact mass of interest)"
        leveltext = "{:s} ({:s})".format(level, level_exp)

        #inventory.append("L4 (Unequivocal precursor ion formula)")
        inventory.append(leveltext)
        inventory.append(hitformula)
        inventory.append(ionform)
        #inventory.append("diff (ppm) {:.4f}".format(diff))

        inventory.append(level_exp2.format(diff))

        ei = ddatoolbox.paths.getsettings("RESOURCES", "ei_fm_def")
        ei_exp = ddatoolbox.paths.getsettings("RESOURCES", "ei_fm_exp")
        

        evidencelist = "[{:s}]".format(",".join([ei]))
        #inventory2.append("[FM]") 
        inventory2.append(evidencelist) 

        hitrecord["formula"] = hitformula
        hitrecord["ionform"] = ionform
        rankinglist = []
        for method, temp in anno["tg"][hitformula].items():
            if method in usedmethods:
                texttemp = "{} (rank {:d})".format(method, int(temp["rank"]))
                rankinglist.append(texttemp)

        #proposetext = "FM (Formula/mass evidence) q-Value-like = {:.6f}, ".format(score)
        proposetext = "{} ({}) {} = {:.6f}, ".format(ei, ei_exp, qlikevalue, score)

        proposetext = proposetext + ", ".join(rankinglist)
        inventory2.append(proposetext) 

    #print("; ".join(inventory))
    annotationtext = "; ".join(inventory)
    annotationtext2 = "; ".join(inventory2)
    return(annotationtext, hitrecord, annotationtext2)

def hitsiriusconf(anno, method = "SIRIUS6_Confidence:biodatabase", thres = 0.64):
    """Filters SIRIUS candidates based on a confidence score threshold.

    Note:
        This helper is not currently called elsewhere in this module; it is
        kept as a standalone utility for ad-hoc SIRIUS-confidence filtering.

    Args:
        anno (dict): Annotation data parsed by `MassBankRecord.get_proposed`.
        method (str, optional): The SIRIUS method key. Defaults to "SIRIUS6_Confidence:biodatabase".
        thres (float, optional): The confidence score threshold. Defaults to 0.64.

    Returns:
        list: A list of InChIKeys that pass the threshold.
    """
    #
    # SIRIUS6 Confidence > 0.64, top10 hits
    #
    SIRIUS_hits = []
    if method in anno["methodtg"]:
        for inchikey, datatemp in anno["methodtg"][method].items():
            score = datatemp["score"]
            if score > thres:
                SIRIUS_hits.append(inchikey)
    return(SIRIUS_hits)

def checklevel3support(anno, inchikey, InChIkeydb, supportmethods = ['MetFrag:MS-FINDER_DB','MS-FINDER:MS-FINDER_DB','CFM-ID:MS-FINDER_DB']):
    """Evaluates supporting evidence for a Level 3 candidate from multiple tools and databases.

    Args:
        anno (dict): Annotation data parsed by `MassBankRecord.get_proposed`.
        inchikey (str): The InChIKey of the candidate to check.
        InChIkeydb (dict): Local database of compound metadata.
        supportmethods (list, optional): List of methods to consider as support.
            Defaults to ['MetFrag:MS-FINDER_DB', 'MS-FINDER:MS-FINDER_DB', 'CFM-ID:MS-FINDER_DB'].

    Returns:
        dict: A dictionary with the following keys:
            - "ranking" (list): (method, score, rank) tuples for candidates
              ranked <= 10 by a `supportmethods` tool.
            - "entropy" (list): (method, score, rank) tuples for GNPS
              entropy-similarity scores > 0.7.
            - "similarity" (list): (method, score, rank) tuples for GNPS
              dot-product similarity scores > 0.7.
            - "hmbd" (str): HMDB ID looked up from `InChIkeydb`, or "N/A".
            - "rankingtext" (str): Human-readable summary of "ranking".
            - "entropytext" (str): Human-readable "SL (Spectral library
              evidence) ..." summary built from "entropy".
            - "similaritytext" (str): Human-readable "SL (Spectral library
              evidence) ..." summary built from "similarity".
            - "hmbdtext" (str): Human-readable "CX (Contextual evidence) ..."
              summary built from "hmbd".
    """
    #supportmethods = ['MetFrag:MS-FINDER_DB','MS-FINDER:MS-FINDER_DB','CFM-ID:MS-FINDER_DB']
    rankinglist = []

    entropymethods = ['Entropy:GNPS']
    entropylist = []

    similaritymethods = ['SIMILARITY:GNPS']
    similaritylist = []



    for method, temp in anno["tg"][inchikey].items():
        if method in supportmethods:
            rank = temp["rank"]
            score = temp["score"]
            if rank <= 10:
                rankinglist.append((method, score, rank))
        if method in similaritymethods:
            rank = temp["rank"]
            score = temp["score"]
            if score >0.7:
                similaritylist.append((method, score, rank))
        if method in entropymethods:
            rank = temp["rank"]
            score = temp["score"]
            if score >0.7:
                entropylist.append((method, score, rank))
    if inchikey in InChIkeydb:
        hmbdid = InChIkeydb[inchikey]['HMDB']
    else:
        hmbdid = "N/A"

    output = {}
    output["ranking"] = rankinglist
    output["entropy"] = entropylist
    output["similarity"] =similaritylist
    output["hmbd"] = hmbdid

    output["rankingtext"] = ["{} (rank {:d})".format(method, rank) for (method, score, rank) in rankinglist]
    output["entropytext"] = ["{} (score {:.3f})".format(method, score) for (method, score, rank) in entropylist]
    output["similaritytext"] = ["{} (score {:.3f})".format(method, score) for (method, score, rank) in similaritylist]
    output["hmbdtext"] = ["{} ({})".format("HMDB", hmdbid) for hmdbid in [hmbdid]]

    output["rankingtext"] = ", ".join(output["rankingtext"])

    ei = ddatoolbox.paths.getsettings("RESOURCES", "ei_spectral_def")
    ei_exp = ddatoolbox.paths.getsettings("RESOURCES", "ei_spectral_exp")
    eitext = "{} ({}) ".format(ei, ei_exp) # "SL (Spectral library) " 
    
    output["similaritytext"] = eitext + ", ".join(output["similaritytext"])
    output["entropytext"] = eitext + ", ".join(output["entropytext"])

    ei = ddatoolbox.paths.getsettings("RESOURCES", "ei_contextual_def")
    ei_exp = ddatoolbox.paths.getsettings("RESOURCES", "ei_contextual_exp")
    eitext = "{} ({}) ".format(ei, ei_exp)
    #output["hmbdtext"] = "CX (Contextual evidence) " + ", ".join(output["hmbdtext"])

    output["hmbdtext"] = eitext + ", ".join(output["hmbdtext"])

    return(output)




def checklevel3(mb2, InChIkeydb, method = "qvalue"):
    """Checks if the record meets the criteria for Level 3-SP confidence (Tentative candidate, structure-proposed).

    The top-ranked compound candidate (by Bayesian q-value-like score)
    qualifies for Level 3-SP when its score is below 0.01. The sub-level
    label is further refined depending on which additional evidence is
    available: in-silico ranking support (IS), GNPS spectral-library
    similarity/entropy (SL), and/or an HMDB biological-context hit (CX).

    Args:
        mb2 (MassBankRecord): The MassBank record object.
        InChIkeydb (dict): Local database of compound metadata (used to
            look up HMDB IDs), keyed by (short) InChIKey.
        method (str, optional): Reserved for future use (not referenced
            inside this function; the q-value-like score is always used).
            Defaults to "qvalue".

    Returns:
        tuple: A tuple containing:
            - annotationtext (str): Formatted confidence string for Level 3-SP
              (empty string if the record does not qualify).
            - hitrecord (dict): Dictionary with keys "inchikey", "hmdb", and
              "level" (e.g. "L3-SP[IS, SL, CX]") describing the winning
              candidate.
            - annotationtext2 (str): Formatted evidence-inventory string
              listing the supporting evidence categories (IS/SL/CX).
    """
    #
    # check level 3
    #
    #COMMENT: CONFIDENCE Tentative candidate (Level 3); SMB Level 3a; Proposed by SIRIUS5_CSI (FDR < 0.05);
    #Similarity by SIMILARITY:GNPS (Dot product > 0.7, 70 spectra); Supported by CFM-ID:MS-FINDER_DB (Ranking 1), MS-FINDER:MS-FINDER_DB (Ranking 1), MetFrag:MS-FINDER_DB (Ranking 5),
    #MetFrag:PubChemLite (Ranking 1); Biological context by HMDB (HMDB0000929)
    #anno = getannodata(mb2, target = "compound")

    qvaluethreshold = 0.05
    anno = mb2.get_proposed(target = "compound")
    
    annotationtext = ""
    annotationtext2 = ""
    inventory = []
    inventory2 = []

    hitrecord = {"inchikey": "", "hmdb": "", "level": ""}

    if not 1 in anno["rank"]:
        return("", hitrecord, "")

    qlikevalue = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue")
    zbramethod = ddatoolbox.paths.getsettings("RESOURCES", "zbramethod")    
    qlikevalue2 = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue2")
 

    methodtext = "q-value_like:Bayesian"
    methodtext = "{:s}:{:s}".format(qlikevalue, zbramethod)
    score1 = 10
    score2 = 10
    score_sirius = 0

    if methodtext in anno["methodrank"]:
        #print(anno["methodrank"][methodtext])
        inchikey1 = anno["methodrank"][methodtext][1]["tg"]
        score1 = anno["methodrank"][methodtext][1]["score"]

    methodtext = "{:s}:{:s}".format(qlikevalue2, zbramethod)
    if methodtext in anno["methodrank"]:
        inchikey2 = anno["methodrank"][methodtext][1]["tg"]
        score2 = anno["methodrank"][methodtext][1]["score"]

    methodtext = "SIRIUS6_Confidence:biodatabase"
    if methodtext in anno["methodrank"]:
        inchikey_sirius = anno["methodrank"][methodtext][1]["tg"]
        score_sirius = anno["methodrank"][methodtext][1]["score"]

    passthres = False


    if min(score1, score2) < qvaluethreshold:
        passthres = True
        if score1 < score2:
            score = score1
            inchikey = inchikey1
        else:
            score = score2
            inchikey = inchikey2
    #elif score_sirius > 0.85:
    elif score_sirius > 1.00:
    

        passthres = True
        inchikey = inchikey_sirius

        score = score_sirius

    if passthres:
        sup = checklevel3support(anno, inchikey, InChIkeydb, supportmethods = ['SIRIUS6_CSI:biodatabase', 'MetFrag:MS-FINDER_DB','MS-FINDER:MS-FINDER_DB','MetFrag:PubChemLite'])

        # Get texts
        level = ddatoolbox.paths.getsettings("RESOURCES", "level3sp_def")
        level_exp = ddatoolbox.paths.getsettings("RESOURCES", "level3sp_exp")
        leveltext = "{:s} ({:s})".format(level, level_exp)
        inventory.append(leveltext)


        ei = ddatoolbox.paths.getsettings("RESOURCES", "ei_insilico_def")
        ei_exp = ddatoolbox.paths.getsettings("RESOURCES", "ei_insilico_exp")
        qlikevalue = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue")

        #proposetext = "IS (In silico evidence) q-Value-like = {:.6f}, ".format(score)
        proposetext = "{} ({}) {} = {:.6f}, ".format(ei, ei_exp, qlikevalue, score)
        proposetext = proposetext+sup["rankingtext"]

        hitrecord["inchikey"] = inchikey

        ei_insilico = ddatoolbox.paths.getsettings("RESOURCES", "ei_insilico_def")
        ei_spectral = ddatoolbox.paths.getsettings("RESOURCES", "ei_spectral_def")
        ei_contextual = ddatoolbox.paths.getsettings("RESOURCES", "ei_contextual_def")

        #if score < qvaluethreshold and len(sup['entropy'])>= 1 and sup['hmbd'] != "N/A":
        if len(sup['entropy'])>= 1 and sup['hmbd'] != "N/A":    

            # Get texts
            inventory = []
            level = ddatoolbox.paths.getsettings("RESOURCES", "level2_def")
            level_exp = ddatoolbox.paths.getsettings("RESOURCES", "level2_exp")
            leveltext = "{:s} ({:s})".format(level, level_exp)
            inventory.append(leveltext)


            inventory.append(inchikey)
            evidencelist = "[{:s}]".format(",".join([ei_insilico, ei_spectral, ei_contextual]))
            inventory2.append(evidencelist)
            #inventory2.append("[IS, SL, CX]")

            inventory2.append(proposetext)
            inventory2.append(sup["entropytext"])
            inventory2.append(sup["hmbdtext"])
            #hitrecord["level"] = "L3-SP[IS, SL, CX]"
            hitrecord["level"] = "{:s}{:s}".format(level, evidencelist)

        #elif score < qvaluethreshold and len(sup['entropy'])>= 1:
        elif len(sup['entropy'])>= 1:    
            # Get texts
            inventory = []
            level = ddatoolbox.paths.getsettings("RESOURCES", "level2_def")
            level_exp = ddatoolbox.paths.getsettings("RESOURCES", "level2_exp")
            leveltext = "{:s} ({:s})".format(level, level_exp)
            inventory.append(leveltext)

            inventory.append(inchikey)
            evidencelist = "[{:s}]".format(",".join([ei_insilico, ei_spectral]))
            inventory2.append(evidencelist)
            #inventory2.append("[IS, SL]")
            inventory2.append(proposetext)
            inventory2.append(sup["entropytext"])
            #hitrecord["level"] = "L3-SP[IS, SL]"
            hitrecord["level"] = "{:s}{:s}".format(level, evidencelist)

        #elif score < qvaluethreshold and sup['hmbd'] != "N/A":
        elif sup['hmbd'] != "N/A":
            inventory.append(inchikey)
            evidencelist = "[{:s}]".format(",".join([ei_insilico, ei_contextual]))
            inventory2.append(evidencelist)
            #inventory2.append("[IS, CX]")
            inventory2.append(proposetext)
            inventory.append(sup["hmbdtext"])
            #hitrecord["level"] = "L3-SP[IS, CX]"
            hitrecord["level"] = "{:s}{:s}".format(level, evidencelist)

        #elif score < qvaluethreshold:
        else:
        
            inventory.append(inchikey)
            evidencelist = "[{:s}]".format(",".join([ei_insilico]))
            inventory2.append(evidencelist)
            #inventory2.append("[IS]")
            inventory2.append(proposetext)
            #hitrecord["level"] = "L3-SP[IS]"
            hitrecord["level"] = "{:s}{:s}".format(level, evidencelist)
    else:
        inventory = []
        inventory2 = []

    annotationtext = "; ".join(inventory)
    annotationtext2 = "; ".join(inventory2)
    return(annotationtext, hitrecord, annotationtext2)

def checklevel3cl(mb2, InChIkeydb, method = "qvalue"):
    """Checks if the record meets the criteria for Level 3-CL confidence (Tentative candidate, class-supported).

    The top-ranked ClassyFire compound-class candidate (by Bayesian
    q-value-like score) qualifies for Level 3-CL when its score is below
    0.01. This level is weaker than Level 3-SP: only the compound class is
    proposed, not a specific structure.

    Args:
        mb2 (MassBankRecord): The MassBank record object.
        InChIkeydb (dict): Local database of compound metadata. Unused in
            this function; kept for a consistent signature with `checklevel3`.
        method (str, optional): Reserved for future use (not referenced
            inside this function; the q-value-like score is always used).
            Defaults to "qvalue".

    Returns:
        tuple: A tuple containing:
            - annotationtext (str): Formatted confidence string for Level 3-CL
              (empty string if the record does not qualify).
            - hitrecord (dict): Dictionary with keys "classyfire" (the
              predicted ClassyFire class) and "level" (e.g. "L3-CL").
            - annotationtext2 (str): Formatted evidence-inventory string
              ("[IS] IS (In silico evidence) ...").
    """

    anno = mb2.get_proposed(target = "class")    
    annotationtext = ""
    annotationtext2 = ""
    inventory = []
    inventory2 = []

    hitrecord = {"classyfire": "", "level": ""}

    if not 1 in anno["rank"]:
        return("", hitrecord, "")

    qlikevalue = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue")
    zbramethod = ddatoolbox.paths.getsettings("RESOURCES", "zbramethod")    


    #methodtext = "q-value_like:Bayesian"
    methodtext = "{:s}:{:s}".format(qlikevalue, zbramethod)
    if methodtext in anno["methodtg"]:
        classyfire = list(anno["methodtg"][methodtext].keys())
        classyfire = classyfire[0]
        score = anno["methodtg"][methodtext][classyfire]["score"]

        # Get texts
        level = ddatoolbox.paths.getsettings("RESOURCES", "level3cp_def")
        level_exp = ddatoolbox.paths.getsettings("RESOURCES", "level3cp_exp")
        #inventory.append("L3-CL (Tentative candidate, Class-supported)")
        leveltext = "{:s} ({:s})".format(level, level_exp)
        inventory.append(leveltext)
        inventory.append(classyfire)
        hitrecord["classyfire"] = classyfire

        ei_insilico = ddatoolbox.paths.getsettings("RESOURCES", "ei_insilico_def")
        ei_exp = ddatoolbox.paths.getsettings("RESOURCES", "ei_insilico_exp")
        qlikevalue = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue")
        if score < 0.05:
            evidencelist = "[{:s}]".format(",".join([ei_insilico]))
            #inventory2.append("[IS]")   
            inventory2.append(evidencelist)          
            #proposetext = "IS (In silico evidence) q-Value-like = {:.6f}".format(score)
            proposetext = "{} ({}) {} = {:.6f}, ".format(ei_insilico, ei_exp, qlikevalue, score)
            inventory2.append(proposetext)
            #hitrecord["level"] = "L3-CL"
            hitrecord["level"] = "{:s}{:s}".format(level, evidencelist)
        else:
            inventory = []
            inventory2 = []

    

    annotationtext = "; ".join(inventory)
    annotationtext2 = "; ".join(inventory2)
    return(annotationtext, hitrecord, annotationtext2)

def curation(datasets, spec_option = "", keyword=""):
    """Main workflow to curate MassBank records by applying multi-level annotation criteria.

    For every averaged MS2 record found in each dataset's step3 output, this
    function evaluates the Level 3-SP, Level 3-CL, Level 4, and Level 5
    confidence criteria (in that priority order, see `checklevel3`,
    `checklevel3cl`, `checklevel4`, `checklevel5`) and applies the
    corresponding curation (compound name/formula/links, curated record
    title, fragment-ion annotation, and a CONFIDENCE comment) to the
    MassBank record. Curated records are written to each dataset's step4
    output folder, and per-dataset counters/logs are written via
    `ddatoolbox.ddatoolbox.writeLog`. The PubChem lookup cache is persisted
    at the end of each dataset via `ddatoolbox.resources.save_pubchemchecked`.

    Args:
        datasets (list): List of dataset names to process. If empty, every
            dataset returned by `ddatoolbox.paths.get_alldataset_names` is
            processed.
        spec_option (str, optional): Option for spectral modification. When
            set to "hitonly", only fragment ions with a matched formula are
            kept in the output spectrum. Defaults to "".
        keyword (str, optional): Filter to process only MS2 filenames
            containing this substring. Defaults to "" (process all files).

    Returns:
        None
    """
    #
    # Reading MsfinderStructureDB to obtain HMDB ids
    #
    softwarepath = "../software/"
    softwareroot = Path(softwarepath).resolve()
    filenametemp = softwareroot/ "MSFINDER"/"Resources"/"MsfinderStructureDB-VS15.esd"
    #print(filenametemp)
    csv_db = ddatoolbox.ddatoolbox.readsmalltable(filenametemp, keywords = ["Title", "Exact mass", "ChEBI","HMDB","KNApSAcK","LipidMAPS", "InChIKeyshort", "InChIkey"], filetype = "tab")
    csv_db = sorted(csv_db, key = lambda x: x["Exact mass"])
    InChIkeydb = {}
    for tempdata in csv_db:
        InChIkey = tempdata["InChIkey"]
        InChIkey_first = InChIkey.split("-")[0]

        InChIkeydb[InChIkey_first] = tempdata
    """
    with open("./res/InChIkeydb.pkl","wb") as f:
        pickle.dump(InChIkeydb, f)


    InChIkeydb = {}
    with open("./res/InChIkeydb.pkl", 'rb') as f:
        InChIkeydb= pickle.load(f)
    """
    compounds, unrealistic = ddatoolbox.resources.load_compoundlist()
    #
    # reading pubchem_checked
    pubchem_checked = ddatoolbox.resources.load_pubchemchecked()
    #
    # for version control
    #
    execdate = datetime.datetime.now().strftime('%Y.%m.%d')
    #
    # for each dataset
    #
    if len(datasets) == 0:
        datasettemp = ddatoolbox.paths.get_alldataset_names()
    else:
        datasettemp = datasets

    print("Curating for", datasettemp, "keyword =", keyword)

    for dataset in datasettemp:
        datapath = ddatoolbox.paths.generate_paths_data(dataset)
        outputpath = datapath["step4path"]+"/"
        accession_stem_ms2 = datapath["accession_stem_ms2"]

        dt_now = datetime.datetime.now().strftime('%Y.%m.%d.%H.%M.%S')
        logfilename =  datapath["logcuration"]
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now], mode = "a")

        filelist = datapath["step3inputlist"]

        recordtemplate = ddatoolbox.massbankio.read_RecordTemplate_dataset(dataset)


        counter = {"level 2" : 0, "level 3" : 0, "level 3cl" : 0,  "level 4" : 0, "level 5" : 0, "total" :0}

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
            # Retrieve basic record information and evaluate every confidence level
            #
            accession = mb2.get_accession() #id

            level5text = checklevel5(mb2)
            method = "qvalue"
            level4text, level4hitrecord, level4inventry = checklevel4(mb2, method = method)
            level3text, level3hitrecord, level3inventry = checklevel3(mb2, InChIkeydb, method = method)
            level3cltext, level3clhitrecord, level3clinventry = checklevel3cl(mb2, InChIkeydb, method = method)
            #
            # Skip automatic L3-SP annotation when the candidate InChIKey is
            # already present in the curated reference compound list
            # (ddatoolbox/compoundlist.txt); such compounds are curated
            # separately and should not be auto-annotated here.
            #
            if level3hitrecord["inchikey"] in compounds and len(level3text) > 0:
                #print("Hit", accession, level3hitrecord["inchikey"])
                level3text = ""
            #
            # Updating Massbank record by record template
            #
            for (fld1, fld2), texttemp in recordtemplate["all"].items():
                if "SP$" in fld1 or "AC$" in fld1 or "CH$" in fld1:
                    if fld1 == "CH$NAME:":
                        continue
                    if fld2 == "SOLVENT":

                        mb2.record[(fld1,fld2)]["data"] = recordtemplate["SOLVENT"]
                        continue
                    if fld2 == "*":
                        mb2.add(fld1, texttemp)
                    else:
                        if fld2 == "SOLVENT":
                            print(fld1, fld2, texttemp)    
                        mb2.add2(fld1, fld2, texttemp)


            
            mb2.add2("COMMENT:", "DATASET_AUTHORS", recordtemplate["dataset_authors"])
            mb2.add2("COMMENT:", "DATASET_TITLE", recordtemplate["dataset_title"])
            mb2.add2("COMMENT:", "DATASET_PUBLICATION", recordtemplate["dataset_publication"])
            mb2.add2("COMMENT:", "DATASET_PUBLICATION_PMID", recordtemplate["dataset_publication_pmid"])
            mb2.add2("COMMENT:", "DATASET_ID", recordtemplate["dataset_id"])

            #mb2.add("SP$SAMPLE:", recordtemplate["sample"])
            #mb2.add2("AC$CHROMATOGRAPHY:", "COLUMN_NAME", recordtemplate["column_name"])
            #
            # Modification of Massbank record
            #
            mb2.clear_anything("COMMENT:", "CONFIDENCE")
            explanationcomment = ""
            warningtext = ""
            #
            # L3-SP
            #
            if len(level3text) > 0:
                level2 = ddatoolbox.paths.getsettings("RESOURCES", "level2_def")
                if level2 in level3text:
                    commenttext =ddatoolbox.paths.getsettings("RESOURCES", "level2_comment")
                    leveltext = level2
                    counter["level 2"] += 1
                else:
                    commenttext = ddatoolbox.paths.getsettings("RESOURCES", "level3sp_comment")
                    leveltext = ddatoolbox.paths.getsettings("RESOURCES", "level3sp_def")
                    counter["level 3"] += 1

                level = level3hitrecord["level"]
                explanationcomment = commenttext.format(level)

                mb2 = modifyclass3(mb2, level3hitrecord, leveltext, pubchem_checked, InChIkeydb, logfilename, unrealistic,  spec_option = spec_option)
                mb2.add2("COMMENT:", "CONFIDENCE", level3text)
                mb2.add2("COMMENT:", "EVIDENCE_INVENTORY", level3inventry)

                if level3hitrecord["inchikey"] in unrealistic:
                    warningtext = "WARNING: This structure is unrealistic due to including {:s}".format(unrealistic[level3hitrecord["inchikey"]])
                else:
                    warningtext = ""


                if len(level3cltext) > 0:
                    name = level3clhitrecord["classyfire"]
                    classes = ["N/A"]
                    classes.append(name)

                    mb2.add("CH$COMPOUND_CLASS:", "; ".join(classes))
                ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, level3text])
                
                counter["total"] += 1

            #
            # L3-CL
            #

            elif len(level3cltext) > 0:
                #
                # If a Level 4 formula candidate is also available, reuse it
                # to annotate fragment-ion formulas even though the overall
                # confidence level assigned to this record is L3-CL.
                #
                if len(level4text) > 0:
                    formula = level4hitrecord["formula"]
                    ionform = level4hitrecord["ionform"]
                    precformula = ddatoolbox.mfsearch.generateprecmf(formula, ionform)


                    annotation, spectra_comp_fragment = fragmentformulasearch(precformula, mb2, logfilename)
                    mb2.add_annotation_from_items(annotation, ["int", "m/z_nocal", "m/z_std(ppm)", "int_std", "cor", "n", "annotation"])

                    #mb2.add("CH$FORMULA:", formula)
                    #exactmass = ddatoolbox.mfsearch.calcmw(formula)
                    #mb2.add("CH$EXACT_MASS:", exactmass)
                    mb2.add2("MS$FOCUSED_ION:","ION_TYPE", ionform)
                else:
                    formula = "N/A"
                    ionform = "N/A"

                level = ddatoolbox.paths.getsettings("RESOURCES", "level3cp_def")
                commenttext = ddatoolbox.paths.getsettings("RESOURCES", "level3cp_comment")
                explanationcomment = commenttext.format(level)

                #explanationcomment = "Automatically annotated at L3-CL."
                #annotationlevel = "Confidence L3-CL"
                annotationlevel = "Confidence {}".format(level) 

                name = level3clhitrecord["classyfire"]
                #classes = ["N/A"]
                #classes.append(name)
                #mb2.add("CH$COMPOUND_CLASS:", "; ".join(classes))

                mb2.add_curatedrecordtitle(name, annotationlevel)

                if spec_option == "hitonly":
                    if len(spectra_comp_fragment) > 0:
                        mb2.add_spectra(spectra_comp_fragment)

                mb2.add2("COMMENT:", "CONFIDENCE", "{}; {}; {}".format(level3cltext, formula, ionform))
                mb2.add2("COMMENT:", "EVIDENCE_INVENTORY", level3clinventry)


                
                ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, level4text])
                #
                # No comment on compound class
                #
                mb2.add("CH$COMPOUND_CLASS:", "N/A")

                counter["level 3cl"] += 1
                counter["total"] += 1


            #
            # L4
            #
            elif len(level4text) > 0:
                formula = level4hitrecord["formula"]
                ionform = level4hitrecord["ionform"]
                precformula = ddatoolbox.mfsearch.generateprecmf(formula, ionform)

                annotation, spectra_comp_fragment = fragmentformulasearch(precformula, mb2, logfilename)
                mb2.add_annotation_from_items(annotation, ["int", "m/z_nocal", "m/z_std(ppm)", "int_std", "cor", "n", "annotation"])
                name = "{} {}".format(formula,ionform)

                level = ddatoolbox.paths.getsettings("RESOURCES", "level4_def")
                commenttext = ddatoolbox.paths.getsettings("RESOURCES", "level4_comment")
                explanationcomment = commenttext.format(level)

                #explanationcomment = "Automatically annotated at L4."
                #annotationlevel = "Confidence L4"
                annotationlevel = "Confidence {}".format(level) 


                #mb2.add("CH$FORMULA:", formula)
                #exactmass = ddatoolbox.mfsearch.calcmw(formula)
                #mb2.add("CH$EXACT_MASS:", exactmass)
                mb2.add2("MS$FOCUSED_ION:","ION_TYPE", ionform)

                mb2.add_curatedrecordtitle(name, annotationlevel)

                if spec_option == "hitonly":
                    if len(spectra_comp_fragment) > 0:
                        mb2.add_spectra(spectra_comp_fragment)

                mb2.add2("COMMENT:", "CONFIDENCE", level4text)
                mb2.add2("COMMENT:", "EVIDENCE_INVENTORY", level4inventry)

                ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, level4text])
                #
                # No comment on compound class
                #
                mb2.add("CH$COMPOUND_CLASS:", "N/A")

                counter["level 4"] += 1
                counter["total"] += 1


            #
            # L5
            #
            elif len(level5text) > 0:
                mb2.add2("COMMENT:", "CONFIDENCE", level5text)

                level = ddatoolbox.paths.getsettings("RESOURCES", "level5_def")
                commenttext = ddatoolbox.paths.getsettings("RESOURCES", "level5_comment")
                explanationcomment = commenttext.format(level)                
                #explanationcomment = "Automatically annotated at L5."

                ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, level5text])
                annotationlevel = "Confidence {}".format(level) 
                precmz, sd = mb2.get_precmz()
                name = "m/z {}".format(precmz)

                mb2.add_curatedrecordtitle(name, annotationlevel)


                counter["level 5"] += 1
                counter["total"] += 1
                #print(accession, "Level 5", counter["level 5"],counter["total"])


                #
                # No comment on compound class
                #
                mb2.add("CH$COMPOUND_CLASS:", "N/A")

            else:
                mb2.add2("COMMENT:", "CONFIDENCE", "Not available")
                explanationcomment = "Automatically annotated at no level."
                ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, "not available"])

                #
                # No comment on compound class
                #
                mb2.add("CH$COMPOUND_CLASS:", "N/A")
                
            logtext = "Level 2/3sp/3cp/4/5/total {:>6d}/{:>6d}/{:>6d}/{:>6d}/{:>6d}/{:>6d}".format(counter["level 2"], counter["level 3"], counter["level 3cl"], counter["level 4"], counter["level 5"],counter["total"])
            print('\r', "{} {:<6s}; {}".format(dataset, accession, logtext), end = "")
            #
            # add accuracy
            #
            ms1thres = mb2.get_anything("COMMENT:","MS1_M/Z_AVERAGED_THRES")
            ms2thres = mb2.get_anything("COMMENT:","MS2_M/Z_AVERAGED_THRES")

            mb2.add2("AC$MASS_SPECTROMETRY:", "MASS_ACCURACY", "MS1 {} ppm, MS2 {} ppm".format(ms1thres, ms2thres))

            #
            # comment update
            #
            dataset_file_num = int(mb2.get_anything("COMMENT:", "DATASET_FILE_NUM"))
            num_averaged_spec = int(mb2.get_anything("COMMENT:", "NUM_AVERAGED_SPEC"))
            comment1 = "This is an averaged mass spectra data created by averaging {:d} product ion spectra among {:d} DDA data files in the dataset described below."
            comment1 = comment1.format(num_averaged_spec, dataset_file_num)





            mb2.clear_anything("COMMENT:", "*")
            mb2.add("COMMENT:", explanationcomment)
            mb2.add("COMMENT:", comment1)
            if len(warningtext) >= 2:
                mb2.add("COMMENT:", warningtext)
            #
            # Record curation date (for version tracking)
            #
            mb2.add_date(date = execdate, target = "Curated")
            #
            # Cleaning up Formula
            #
            formula = mb2.get_formula()
            formula = re.sub('[-+][0-9]*', "", formula)
            mb2.add("CH$FORMULA:", formula)
            #
            # Add project and author information
            #
            projecttext = ddatoolbox.paths.getsettings("RESOURCES", "project")
            authorstext = ddatoolbox.paths.getsettings("RESOURCES", "authors")
            publicationtext = ddatoolbox.paths.getsettings("RESOURCES", "publication")

            mb2.add("PROJECT:", projecttext)
            mb2.add("AUTHORS:", authorstext)
            mb2.add("PUBLICATION:", publicationtext)
            #
            # Write out the curated record file
            #
 

            mb2.show(path = outputpath)
        ddatoolbox.resources.save_pubchemchecked(pubchem_checked)

        dt_now = datetime.datetime.now()
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Finished", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S')])
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Level 3", accession_stem_ms2, counter["level 3"]])
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Level 4", accession_stem_ms2, counter["level 4"]])
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Level 5", accession_stem_ms2, counter["level 5"]])
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Total", accession_stem_ms2, counter["total"]])
        print("")
        #print(dataset, logtext)
def test():
    """Test function for debugging and verifying specific annotation logic and database retrieval.

    Returns:
        None
    """
    #mb2 = ddatoolbox.massbankio.MassBankRecord()
    #mb2.read_MassBankRecord("../data/HSA009/step4/MSBNK-SMB_Measured-HSA009P0227001.txt")
    #logfilename = "temp.txt"
    mf = "C11H12N2O2"
    ionformpdcompound = "[M+Na]+"

    pdproperty = ddatoolbox.resources.get_pubchemdatafrominchikey("KRHNXNZBLHHEIU-UHFFFAOYSA-N")
    print(pdproperty)


    precmf = ddatoolbox.mfsearch.generateprecmf(mf, ionformpdcompound)
    print(precmf)


    #mb3 = fragmentformulasearch(precmf, mb2, logfilename)


    #ACCESSION: MSBNK-SMB_Measured-HSA009P0227001
    # Level 4
    #MSBNK-SMB_Measured-HSA001P0213002.txt
    # Level 4
    #MSBNK-SMB_Measured-HSA001P0205016.txt



if __name__ == '__main__':
    # Usage: python DDACuration.py <dataset> [keyword]
    # Falls back to a hard-coded default (HSA001, keyword = "") when
    # called with no arguments.
    args = sys.argv
    if len(args) == 3:
        dataset = args[1]
        keyword = int(args[2])
        datasets = [dataset]
    else:
        dataset = ["HSA001"]
        keyword = "HSA001P0205016"

    #datasets = ["HSA045"]
    #level 5 HSA001P0205017
    #level 4 HSA001P0205009
    #level 3cl HSA001P0214003
    #level 3sp  HSA001P0205016
    #level 2  HSA001P0205016
    # unrealistic ['HSA045P0491004', 'HSA059P0491000']
    curation(datasets, spec_option = "", keyword = "")





