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

import ddatoolbox
import argparse
import glob, copy, pickle, csv, itertools, sys, datetime
import os, re, shutil


description = """
ShinMassBank pipeline\n
Example 1:\n
python DDApipeline.py HSA001 HSA002 -a -p -anno MS-FINDER MetFragMSFINDER MetFragPubchemLite -r -c -b 20 -t 5\n
\n
Example 2: pw.txt is a single line text file describing password for your sirius account\n
python DDApipeline.py HSA001 HSA002 -anno SIRIUS5 SIRIUS6 < pw.txt\n

"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument("datasets", nargs='*', default=[])    # 必須の引数を追加
parser.add_argument("-v", "--hiveconvert", action='store_true', help='Conversion of raw data files (".raw", ".lcd", ".wiff", ".d", ".mzML", ".mzXML") to hive format. ')    # 必須の引数を追加
parser.add_argument("-a", "--averaging", action='store_true', help='Execute averaging of MS2 data')    # 必須の引数を追加
parser.add_argument("-b", "--batchsize", default=5, help='Batch size of unit mass for averaging')    # 必須の引数を追加
parser.add_argument("-t", "-n", "--thread", default=1, help='Number of CPU threads for averaging')    # 必須の引数を追加


parser.add_argument("-p", "--processing", action='store_true', help='Execute processing of averaged spec for re-calibraion')    # 必須の引数を追加

parser.add_argument("-anno", "--annotation", nargs='*', default=[""], help='Execute Annotation, SIRIUS6 MS-FINDER MetFragMSFINDER MetFragPubChemLite SimilarityGNPS SimilarityMassBank')    # 必須の引数を追加
parser.add_argument("-s", "--reset", action='store_true', help='Reset annotation result')    # 必須の引数を追加

parser.add_argument("-cfmid", "--containerlist", nargs='*', default=[""], help='List of container name for CFM-ID')    # 必須の引数を追加
parser.add_argument("-sirius", "--siriusid", default="", help='SIRIUS accound id')    # 必須の引数を追加
parser.add_argument("-r", "--retrieve", action='store_true', help='Retrieve all annotation results')    # 必須の引数を追加
parser.add_argument("-f", "--fdr", action='store_true', help='FDR Estimation by Zipf-Bayesian Rank Aggregation')    # 必須の引数を追加
parser.add_argument("-c", "--curation", action='store_true', help='Curate annotation level')    # 必須の引数を追加

parser.add_argument("-l", "--finalize", action='store_true', help='Generation of final records and mzTab')    # 必須の引数を追加


parser.add_argument("-k", "--keyword", default="", help='Keyword in accession for curation test purpose')    # 必須の引数を追加


args = parser.parse_args()    # 4. 引数を解析

datasets = args.datasets
batchsize = int(args.batchsize)
numcpu = int(args.thread)
containerlist = args.containerlist
siriusid = args.siriusid
reset = args.reset
keyword = args.keyword

if args.hiveconvert == True:
    import DDAAveraging
    for dataset in datasets:
        DDAAveraging.hiveconvert(dataset)

if args.averaging == True:
    import DDAAveraging
    print(datasets, batchsize, numcpu)
    for dataset in datasets:
        DDAAveraging.averaging(dataset, batchsize=batchsize, numcpu=numcpu)

if args.processing == True:
    import DDAProcessing
    DDAProcessing.processing(datasets)

if not args.annotation == [""]:
    import DDAAnnotation
    for software in args.annotation:
        if software == "CFM-ID":
            print(software, datasets, containerlist, siriusid)
            DDAAnnotation.annotation(datasets, software, containerlist = containerlist, reset = reset)

        elif "SIRIUS" in software:
            print(software, datasets, siriusid)
            DDAAnnotation.annotation(datasets, software, containerlist = containerlist, siriusid = siriusid, reset = reset)
        else:
            print(software, datasets)
            DDAAnnotation.annotation(datasets, software, reset = reset)


if args.retrieve == True:
    import DDARetrieve
    DDARetrieve.retrieve(datasets, keyword=keyword, step = "first")

if args.fdr == True:
    import DDARetrieve
    import DDAFdr
    print(datasets, keyword)

    parameter = DDAFdr.setparameter()
    confidence_threshold = parameter["confidence_threshold"]
    formulaconfidence_threshold = parameter["formulaconfidence_threshold"]

    thresholds = {}
    thresholds["compound"] = parameter["confidence_threshold"]
    thresholds["formula"] = parameter["formulaconfidence_threshold"]
    thresholds["class"] = parameter["class_threshold"]

    rankinglimit = parameter["rankinglimit"]
    rankinglimit_formula = parameter["rankinglimit_formula"]
    rankinglimit_class = parameter["rankinglimit_class"]

    inputstep = parameter["inputstep"]
    reset = parameter["reset"]

    abbreviations = {}
    abbreviations["compound"] = parameter["abbreviations"]
    abbreviations["formula"] = parameter["abbreviations_formula"]
    abbreviations["class"] = parameter["abbreviations_class"]

    combination = parameter["combination"]
    combination_formula = parameter["combination_formula"]
    combination_class = parameter["combination_class"]

    DDAFdr.zipfannotation(datasets, abbreviations["compound"], thresholds["compound"], target = "compound", keyword = keyword, inputstep = inputstep, rankinglimit = rankinglimit, combination = combination, reset = True)
    DDAFdr.zipfannotation(datasets, abbreviations["formula"], thresholds["formula"], target = "formula", keyword = keyword, inputstep = inputstep, rankinglimit = rankinglimit_formula, combination = combination_formula, reset = reset)
    DDAFdr.zipfannotation(datasets, abbreviations["class"], thresholds["class"], target = "class", keyword = keyword, inputstep = inputstep, rankinglimit = rankinglimit_class, combination = combination_class, reset = reset)
    DDARetrieve.retrieve(datasets, keyword=keyword)

if args.curation == True:
    import DDACuration

    DDACuration.curation(datasets, keyword=keyword)

if args.finalize == True:
    import DDAFinalize

    for dataset in datasets:
        #print("Generating records for MassBank Human", dataset)
        #DDAFinalize.finalize(dataset, keyword=keyword)
        print("Generating mzTab-m", dataset)
        DDAFinalize.generate_mztab(dataset, keyword=keyword)