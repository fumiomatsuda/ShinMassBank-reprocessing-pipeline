# -------------------------------------------------------------------------------
# Name:        DDARetrieve
# Purpose:     Collect the annotation results produced by every annotation tool
#              (SIRIUS 5/6, MS-FINDER, CFM-ID, MetFrag, GNPS/MassBank similarity
#              search, and the Bayesian Zipf/FDR pipeline) and merge them, sorted
#              by score, back into each dataset's MassBank record files.
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
# -------------------------------------------------------------------------------

import ddatoolbox
import glob, copy, pickle, csv, itertools, sys, datetime
from pathlib import Path
import os, re


def retrieve(datasets, keyword="", step="notfirst"):
    """Merge every annotation tool's output back into each record's PROPOSED_* comments.

    For each dataset, reads the ``annotation/<software>/output.txt`` result table of
    every configured annotation tool (SIRIUS 5/6, MS-FINDER, CFM-ID, MetFrag, GNPS/
    MassBank similarity search, and, when ``step != "first"``, the Bayesian
    (``Baysian``) and 2nd-rank FDR (``FDR_2ndrank``) post-processing outputs),
    indexed by ``(accession, target_type)``. For every MassBank record under
    ``step2``, existing ``PROPOSED_MF``, ``PROPOSED_MF_FROM_PRECMZ``,
    ``PROPOSED_FROM_SPEC_SIMILARITY``, ``PROPOSED_COMPOUND``, and ``PROPOSED_CLASS``
    comments are cleared and rebuilt: candidates from each tool are sorted by score
    (direction looked up per ``(target_type, method)`` in ``sortingorder``,
    defaulting to ascending) and appended as new comment lines, with duplicate
    entries skipped. Compound-level hits from a fixed set of "class-worthy" methods
    are additionally mapped to ClassyFire compound classes (via
    ``ddatoolbox.resources.load_inchikey2classyfire``) and recorded as
    ``PROPOSED_CLASS`` entries, ranked by first occurrence. Updated records are
    stamped with the current date and written out to the ``step3`` folder.

    Args:
        datasets (list[str]): Dataset identifiers to process.
        keyword (str, optional): Filter string; only record files whose filename
            contains ``keyword`` are processed. If empty, all files are processed.
            Defaults to ``""``.
        step (str, optional): When ``"first"``, only the primary annotation tools
            (SIRIUS 5/6, MS-FINDER, CFM-ID, MetFrag, similarity search) are read.
            Any other value additionally includes the ``Baysian`` and
            ``FDR_2ndrank`` post-processing outputs. Defaults to ``"notfirst"``.

    Returns:
        None
    """



    InChIkeydb = ddatoolbox.resources.load_inchikey2classyfire()

    execdate = datetime.datetime.now().strftime("%Y.%m.%d")

    plikevalue = ddatoolbox.paths.getsettings("RESOURCES", "plikevalue")
    qlikevalue = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue")
    zbramethod = ddatoolbox.paths.getsettings("RESOURCES", "zbramethod")   


    if step == "first":
        annotation_software = [
            "SIRIUS5",
            "SIRIUS6",
            "MS-FINDER",
            "CFM-ID",
            "MetFragMSFINDER",
            "MetFragPubChemLite",
            "SimilarityGNPS",
            "SimilarityMassBank",
        ]
    else:
        annotation_software = [
            "SIRIUS5",
            "SIRIUS6",
            "MS-FINDER",
            "CFM-ID",
            "MetFragMSFINDER",
            "MetFragPubChemLite",
            "SimilarityGNPS",
            "SimilarityMassBank",
            zbramethod,

        ]

    #
    # Sort direction for each software
    #

    sortingorder = {
        ("mf", "SIRIUS5_Zodiac"): "Descending",
        ("mf", "SIRIUS5_Sirius"): "Descending",
        ("mf", "MS-FINDER"): "Descending",
        ("mf", "SIRIUS6_Zodiac"): "Descending",
        ("mf", "SIRIUS6_Sirius"): "Descending",
        ("mf", "Aggregated_ranking"): "Descending",
        ("mf", plikevalue): "Descending",
        ("mf", qlikevalue): "Descending",
        ("compound", "SIRIUS6_Confidence"): "Descending",
        ("compound", "SIRIUS6_CSI"): "Descending",
        ("compound", "SIRIUS5_Confidence"): "Descending",
        ("compound", "SIRIUS5_CSI"): "Descending",
        ("compound", "MS-FINDER"): "Descending",
        ("compound", "CFM-ID"): "Descending",
        ("compound", "MetFrag"): "Descending",
        ("compound", "MetFragMSFINDER"): "Descending",
        ("compound", "MetFragPubChemLite"): "Descending",
        ("compound", "SimilarityGNPS"): "Descending",
        ("compound", "SimilarityMassBank"): "Descending",
        ("compound", "Aggregated_ranking"): "Descending",
        ("compound", "Aggregated_ranking_2"): "Descending",
        ("compound", "Aggregated_ranking_3"): "Descending",
        ("compound", plikevalue): "Descending",
        ("compound", qlikevalue): "Descending",
        ("compound", plikevalue+"_2"): "Descending",
        ("compound", qlikevalue+"_2"): "Descending",
        ("compound", plikevalue+"_3"): "Descending",
        ("compound", qlikevalue+"_3"): "Descending",
        ("class", "SIRIUS6_CSI"): "Descending",
        ("class", "SIRIUS5_CSI"): "Descending",
        ("class", "MS-FINDER"): "Descending",
        ("class", "MetFragMSFINDER"): "Descending",
        ("class", "MetFragPubChemLite"): "Descending",
        ("class", "SimilarityGNPS"): "Descending",
        ("class", "SimilarityMassBank"): "Descending",
        ("class", "Aggregated_ranking"): "Descending",
        ("class", plikevalue): "Descending",
        ("class", qlikevalue): "Descending",
    }

    classtargetmethods = [
        "Entropy:GNPS",
        "MS-FINDER:MS-FINDER_DB",
        "MetFrag:MS-FINDER_DB",
        "MetFrag:PubChemLite",
        "SIRIUS6_CSI:biodatabase",
        "SIRIUS5_CSI:biodatabase",
    ]

    for dataset in datasets:
        paths = ddatoolbox.paths.generate_paths_data(dataset)

        inputpath = paths["step2input"]
        outputpath = paths["step3path"]
        annotationpath = paths["annotation"] + "/"
        logpath = paths["logpath"]


        files = glob.glob(inputpath)
        #
        # Reading output.txt in /annotation/
        #
        annotation_data = {}

        for software in annotation_software:
            annotationfile = annotationpath + software + "/output.txt"
            if not os.path.isfile(annotationfile):
                continue
            with open(annotationfile, encoding="UTF-8") as f:
                reader = csv.reader(f, delimiter="\t")
                for i, row in enumerate(reader):
                    if len(row) < 8:
  
                        continue
                    id = row[0]  # redord id
                    type = row[1]  # compound
                    scoring = row[2]  # MetFrag
                    db = row[3]  # PubChem
                    score = float(row[4])  # 1.0
                    result = row[5]  # PSZFDJIBPJBVSH
                    ionform = row[6]  # [M+H]+
                    name = "_".join(row[7:])
                    name = name.split("\t")[0]  # (2S,3R)-2-aminotridec-4-yne-1,3-diol
                    name = name.replace('"', "")

                    #
                    # For tie score of SIRIUS
                    #
                    if ionform == "[M+H]+" and scoring == "SIRIUS5_Sirius":
                        score = score + 0.0001
                    if ionform == "[M+H]+" and scoring == "SIRIUS6_Sirius":
                        score = score + 0.0001
                    if ionform == "[M+H]+" and scoring == "SIRIUS5_Zodiac":
                        score = score + 0.0001
                    if ionform == "[M+H]+" and scoring == "SIRIUS6_Zodiac":
                        score = score + 0.0001

                    if db == "MsfinderDB":
                        db = "MS-FINDER_DB"
                    if db == "PubchemLite":
                        db = "PubChemLite"

                    scoredb = (scoring, db)
                    idtype = (id, type)

                    if idtype not in annotation_data:
                        annotation_data[idtype] = {}

                    if scoredb not in annotation_data[idtype]:
                        annotation_data[idtype][scoredb] = []
                    temphash = {
                        "score": score,
                        "result": result,
                        "ionform": ionform,
                        "name": name,
                    }
                    annotation_data[idtype][scoredb].append(temphash)

        # List record files under /step2/
        #
        filelist = []
        for file in files:
            if not ".txt" in file:
                continue
            filelist.append(file)
        # Read record files
        #

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
            # Get record info
            #
            accession = mb2.get_accession()  # get the accession ID
            polarity = mb2.get_polarity()  # get the ion polarity
            precmz, sd = mb2.get_precmz()

            # Get remaining fields directly


            mb2.clear_anything("COMMENT:", "PROPOSED_MF_FROM_PRECMZ")
            mb2.clear_anything("COMMENT:", "PROPOSED_FROM_SPEC_SIMILARITY")
            mb2.clear_anything("COMMENT:", "PROPOSED_COMPOUND")

            checkduplicate = {}
            #
            # Molecular formula
            #
            type = "mf"
            if (accession, type) in annotation_data:
                scoredblist = sorted(
                    annotation_data[(accession, type)].keys(), key=lambda x: x
                )
                for scoredb in scoredblist:
                    (scoringmethod, db) = scoredb
                    sortdirection = "Ascending"
                    if (type, scoringmethod) in sortingorder:
                        sortdirection = sortingorder[(type, scoringmethod)]

                    hittedlist = annotation_data[(accession, type)][scoredb]
                    if sortdirection == "Ascending":
                        hittedlistsorted = sorted(hittedlist, key=lambda x: x["score"])
                    else:
                        hittedlistsorted = sorted(
                            hittedlist, key=lambda x: x["score"], reverse=True
                        )

                    for hitted in hittedlistsorted:
                        software = scoringmethod
                        formula = hitted["result"]
                        score = hitted["score"]
                        ionform = hitted["ionform"]
                        name = hitted["name"]

                        precmf = ddatoolbox.mfsearch.generateprecmf(formula, ionform)
                        exactprecmz = ddatoolbox.mfsearch.calcmw(
                            precmf, polarity=polarity
                        )
                        delta = (precmz - exactprecmz) / exactprecmz * 1000000

                        thiskey = ("COMMENT:", "PROPOSED_MF")
                        formattext = mb2.record[thiskey]["format"]

                        if software in [
                            "Aggregated_ranking",
                            "Aggregated_ranking_2",
                            "Aggregated_ranking_3",
                            plikevalue,
                            qlikevalue,
                            plikevalue+"_2",
                            qlikevalue+"_2",
                            plikevalue+"_3",
                            qlikevalue+"_3",
                        ]:
                            formattext = formattext.replace("{:.3f}", "{:.7f}")

                        text = formattext.format(
                            software, score, formula, ionform, delta
                        )
                        if (software, score, formula, ionform, delta) in checkduplicate:
                            continue
                        checkduplicate[(software, score, formula, ionform, delta)] = 0

                        mb2.add2("COMMENT:", "PROPOSED_MF", text)

            #
            # Attach compound information
            #
            # Compile the information collected from output.txt into comment text
            type = "compound"
            if (accession, type) in annotation_data:
                scoredblist = sorted(
                    annotation_data[(accession, type)].keys(), key=lambda x: x
                )
                for scoredb in scoredblist:
                    (scoringmethod, db) = scoredb
                    sortdirection = "Descending"
                    if (type, scoringmethod) in sortingorder:
                        sortdirection = sortingorder[(type, scoringmethod)]

                    hittedlist = annotation_data[(accession, type)][scoredb]
                    if sortdirection == "Ascending":
                        hittedlistsorted = sorted(hittedlist, key=lambda x: x["score"])
                    else:
                        hittedlistsorted = sorted(
                            hittedlist, key=lambda x: x["score"], reverse=True
                        )

                    for hitted in hittedlistsorted:
   
                        software = scoringmethod
                        database = db
                        inchikey = hitted["result"]
                        score = float(hitted["score"])
                        ionform = hitted["ionform"]
                        name = hitted["name"]

                        if software == "CFM=ID" and score <= 0.0:
                            continue
                        #
                        # Spectral Match to L-Tryptophan from NIST14::C11H12N2O2::131::0.9897::0.8571::0.9145
                        # Spectral Match to L-Tryptophan from NIST14 C11H12N2O2
                        if software == "SIMILARITY" or software == "Entropy":
                            name = name.split("::")[:2]

                            name = "_".join(name)

                        #if "A-N" in inchikey:
                        #    continue

                        
                        thiskey = ("COMMENT:", "PROPOSED_COMPOUND")
                        formattext = mb2.record[thiskey]["format"]

                        if software in [
                            "Aggregated_ranking",
                            plikevalue,
                            qlikevalue,
                        ]:
                            formattext = formattext.replace("{:.3f}", "{:.7f}")

                        text = formattext.format(
                            software, database, score, inchikey, ionform, name
                        )

                        if (
                            software,
                            database,
                            score,
                            inchikey,
                            ionform,
                            name,
                        ) in checkduplicate:
                            continue
                        checkduplicate[
                            (software, database, score, inchikey, ionform, name)
                        ] = 0

                        mb2.add2("COMMENT:", "PROPOSED_COMPOUND", text)
                    #
                    # class
                    #
                    method = ":".join([scoringmethod, db])

                    if not method in classtargetmethods:
                        continue

                    classrecord = {}
                    counter = set()
                    rank = 0
                    for hitted in hittedlistsorted:
                        rank += 1
                        software = scoringmethod
                        database = db
                        inchikey = hitted["result"]
                        inchikey_short = inchikey.split("-")[0]
                        score = float(hitted["score"])
                        ionform = hitted["ionform"]
                        name = hitted["name"]
                        if software == "CFM=ID" and score <= 0.0:
                            continue
                        if not inchikey_short in InChIkeydb:
                            continue
                        classyfire = InChIkeydb[inchikey_short]
                        if not classyfire in counter:
                            counter.add(classyfire)
                            classrank = len(counter)
                            classrecord[classrank] = {
                                "classyfire": classyfire,
                                "score": score,
                                "rank": rank,
                                "inchikey": inchikey_short,
                                "ionform": ionform,
                                "name": name,
                            }

                    for rank in classrecord.keys():
                        classyfire = classrecord[rank]["classyfire"]
                        thisrank = classrecord[rank]["rank"]
                        inchikey = classrecord[rank]["inchikey"]
                        name = classrecord[rank]["name"]
                        ionform = classrecord[rank]["ionform"]
                        score = classrecord[rank]["score"]
                        # outputtemp = [accession, "class", method, "ClassyFire", rank, classyfire, inchikey, thisrank]
                        # outputall.append(outputtemp)
                        thiskey = ("COMMENT:", "PROPOSED_CLASS")
                        formattext = mb2.record[thiskey]["format"]
                        text = formattext.format(
                            software,
                            database,
                            score,
                            classyfire,
                            ionform,
                            "_".join([inchikey, name]),
                        )

                        if text in checkduplicate:
                            continue
                        checkduplicate[text] = 0

                        mb2.add2("COMMENT:", "PROPOSED_CLASS", text)

            #
            # Attach class information
            #
            # Compile the information collected from output.txt into comment text
            type = "class"
            if (accession, type) in annotation_data:
                scoredblist = sorted(
                    annotation_data[(accession, type)].keys(), key=lambda x: x
                )
                for scoredb in scoredblist:
                    (scoringmethod, db) = scoredb
                    
                    if scoringmethod not in ["Aggregated_ranking", plikevalue, qlikevalue]:
                        continue
                    
 
                    sortdirection = "Descending"
                    if (type, scoringmethod) in sortingorder:
                        sortdirection = sortingorder[(type, scoringmethod)]

                    hittedlist = annotation_data[(accession, type)][scoredb]
                    if sortdirection == "Ascending":
                        hittedlistsorted = sorted(hittedlist, key=lambda x: x["score"])
                    else:
                        hittedlistsorted = sorted(
                            hittedlist, key=lambda x: x["score"], reverse=True
                        )

                    for hitted in hittedlistsorted:
                        software = scoringmethod
                        database = db
                        inchikey = hitted["result"]
                        score = float(hitted["score"])
                        ionform = hitted["ionform"]
                        name = hitted["name"]

                        if "A-N" in inchikey:
                            continue

                        thiskey = ("COMMENT:", "PROPOSED_CLASS")
                        formattext = mb2.record[thiskey]["format"]

                        if software in [
                            "Aggregated_ranking",
                            plikevalue,
                            qlikevalue,
                        ]:
                            formattext = formattext.replace("{:.3f}", "{:.7f}")

                        text = formattext.format(
                            software, database, score, inchikey, ionform, name
                        )

                        if text in checkduplicate:
                            continue
                        checkduplicate[text] = 0

                        mb2.add2("COMMENT:", "PROPOSED_CLASS", text)

            #
            # For version tracking
            #

            mb2.add_date(date=execdate, target="Annotated")
            #
            # Write the file out
            #
            mb2.show(path=outputpath)
            print("\r", "Retrieve:", ms2filename, end="")
    print("")


if __name__ == "__main__":
    # Usage: python DDARetrieve.py <dataset1> [dataset2 ...]
    # Falls back to a hard-coded default (HSA001, filtered by keyword "205016")
    # when called with no arguments.
    args = sys.argv
    if len(args) >= 2:
        datasets = args[1:]
        keyword = ""
    else:
        datasets = ["HSA001old"]
        keyword = "205016"

    retrieve(datasets, keyword=keyword, step="second")
