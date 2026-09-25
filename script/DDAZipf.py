#-------------------------------------------------------------------------------
# Name:        DDAZipf
# Purpose:     Generation of the preestimated Zipf distribution (Bayesian model fitted
#              with PyMC) used by DDAcalcfdr.py for integrated ranking and FDR
#              estimation of MassBank records.
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import ddatoolbox
import glob, copy, pickle, csv, itertools, datetime, argparse
from pathlib import Path
import re

import numpy as np
from matplotlib import pyplot as plt
import pymc as pm
import pandas as pd

from DDAFdr import setparameter, pathsetup



def main(abbreviations, df, rankinglimit, target = "compound"):
    """Build a Bayesian model based on the Zipf distribution and sample the posterior distribution.

    For all pairwise (ordered) combinations of annotation methods, samples the posterior
    distributions of the Zipf-distribution shape parameters (``la``, ``lb``, ``l0``, ``l1``)
    and the database correct-answer rate (``alpha``) from the observed compound-ranking
    frequency data using PyMC (No-U-Turn / NUTS sampler under the hood).

    Args:
        abbreviations (dict): Dictionary mapping annotation-method names to short
            abbreviations used to build PyMC variable names.
            Example: ``{"SIRIUS6_CSI:biodatabase": "S6", ...}``
        df (pandas.DataFrame): DataFrame containing the ranking-count data produced by
            :func:`rankingcount`. Column names follow the pattern
            ``"<firststep>_<nextstep>_all"``, ``"_hit"``, ``"_no"``, ``"_nodb"``.
        rankinglimit (int): Upper limit of ranks (i.e. how many leading rows of ``df``)
            used for fitting the model.

    Returns:
        tuple: A tuple of two elements.

            - trace: PyMC sampling result (``MultiTrace`` object) as returned by
              ``pm.sample()``.
            - tracedict (dict): Dictionary mapping parameter names to lists of sampled
              values. Internal PyMC variables containing ``"interval__"`` (the
              unconstrained reparameterizations of bounded variables) are excluded.
    """
    allobsd = {}
    hitobsd = {}
    noobsd = {}
    nodbobsd = {}
    gtobsd = {}
    methods = list(abbreviations.keys())
    combinations = list(itertools.permutations(methods, 2))

    print("Zipf: Generating posterior distribution for ", methods, "Ranking limit", rankinglimit)

    for firststep, nextstep in combinations:
        cmb = (firststep, nextstep)
        allobsd[cmb], hitobsd[cmb], noobsd[cmb], nodbobsd[cmb],  gtobsd[firststep] = obsdsetup(firststep, nextstep, df)


    with pm.Model() as model:
        X = np.array(range(1,200))
        alpha = pm.Uniform("alpha", lower=0.01, upper = 0.999) # slope (exponent) prior
        sigma = pm.HalfNormal("sigma", sigma=0.01)

        la = {}
        lb = {}
        l0 = {}
        l1 = {}

        zipa = {}
        zipb = {}
        zip0 = {}
        zip1 = {}
        allpred = {}
        hitpred = {}
        nopred = {}
        nodbpred = {}
        gtpred = {}
        Y_hit = {}
        Y_no = {}
        Y_nodb = {}
        Y_all = {}
        Y_gt = {}

        for firststep in methods:
            firststepabr = abbreviations[firststep]
            la[firststep]= pm.Uniform("la_"+firststepabr, lower=-12.0, upper = -1.001) # slope (exponent) prior
            # first hit
            ziptemp = X ** la[firststep]
            zipsum = pm.math.sum(ziptemp)
            zipa[firststep] = ziptemp/zipsum

            gtpred[firststep] = zipa[firststep]

            Y_gt[firststep] = pm.Normal("Y_gt"+firststepabr, mu=gtpred[firststep][:rankinglimit], sigma=sigma, observed=gtobsd[firststep][:rankinglimit])

        for firststep, nextstep in combinations:

            cmb = (firststep, nextstep)
            firststepabr = abbreviations[firststep]
            cmbabr = abbreviations[firststep]+abbreviations[nextstep]

            lb[cmb] = pm.Uniform("lb_"+cmbabr, lower=-12.0, upper = -1.001)
            l0[cmb] = pm.Uniform("l0_"+cmbabr, lower=-8.0, upper = -1.001) # slope (exponent) prior
            l1[cmb] = pm.Uniform("l1_"+cmbabr, lower=-8.0, upper = -1.001) # slope (exponent) prior
            # second hit
            ziptemp = X ** lb[cmb]
            zipsum = pm.math.sum(ziptemp)
            zipb[cmb] = ziptemp/zipsum
            # nohit
            ziptemp = X ** l0[cmb]
            zipsum = pm.math.sum(ziptemp)
            zip0[cmb] = ziptemp/zipsum
            # nodb
            ziptemp = X ** l1[cmb]
            zipsum = pm.math.sum(ziptemp)
            zip1[cmb] = ziptemp/zipsum
            #
            # alpha: Percentage of correct answers in the database
            #    * The percentage of spectral data that contains the necessary information
            # La, zipa: Zipf parameter for first ranking, correct hits
            # Lb, zipb: Zipf parameter for second ranking, correct hit in first ranking
            # L0, zip0: Zipf parameter for second ranking, wrong hit in first ranking, correct answer in database
            # L1, zip1: Zipf parameter for second ranking, wrong hit in first ranking, no correct answer in database
            # allpred: alpha *  hitpred + alpha * nopred + (1 - alpha) * nodbpred
            # hitpred: Zipf_La(1) * Zipf_Lb(rank) First top rank is correct, Second rank n is correct,
            # nopred: (1 - Zipf_La(1)) * (1 - Zipf_Lb(rank)) * Zipf_L0(rank) First top rank is wrong, Second rank n is wrong,
            # nodbpred: Zipf_L1(rank)


            hitpred[cmb] = zipb[cmb]
            nopred[cmb] = (1-zipb[cmb])*zip0[cmb]
            nodbpred[cmb] = zip1[cmb]

            Y_hit[cmb] = pm.Normal("Y_hit"+cmbabr, mu=hitpred[cmb][:rankinglimit], sigma=sigma, observed=hitobsd[cmb][:rankinglimit])
            Y_no[cmb] = pm.Normal("Y_no"+cmbabr, mu=nopred[cmb][:rankinglimit], sigma=sigma, observed=noobsd[cmb][:rankinglimit])
            Y_nodb[cmb] = pm.Normal("Y_nodb"+cmbabr, mu=nodbpred[cmb][:rankinglimit], sigma=sigma, observed=nodbobsd[cmb][:rankinglimit])

            if not target == "formula":
                allpred[cmb] =  alpha * zipa[firststep][0] * zipb[cmb] + alpha * (1-zipa[firststep][0])*(1-zipb[cmb])*zip0[cmb] + (1-alpha) * zip1[cmb]
                Y_all[cmb] = pm.Normal("Y_all"+cmbabr, mu=allpred[cmb][:rankinglimit], sigma=sigma, observed=allobsd[cmb][:rankinglimit])
 

        trace = pm.sample(10000, tune=5000, return_inferencedata=False)

    tracedict = {}
    for varname in trace.varnames:
        if "interval__" in varname:
            continue
        tracedict[varname] = list(trace[varname])

    return(trace, tracedict)

def obsdsetup(firststep, nextstep, df):
    """Extract and normalize the observed rank-count distributions for one method pair.

    For the given ``(firststep, nextstep)`` method combination, reads the four
    rank-indexed count columns from ``df`` and converts each of them into a
    normalized (probability-like) frequency distribution by dividing by the
    appropriate total. These normalized series are the empirical counterparts of
    the Zipf-distribution predictions used in :func:`main`.

    Args:
        firststep (str): Name of the first (upstream) annotation method. Used as the
            first part of the ``"<firststep>_<nextstep>_..."`` column-name prefix.
        nextstep (str): Name of the second (downstream) annotation method. Used as the
            second part of the column-name prefix.
        df (pandas.DataFrame): DataFrame containing the ranking-count data produced by
            :func:`rankingcount`. Must include the columns
            ``"<firststep>_<nextstep>_all"``, ``"_hit"``, ``"_no"``, and ``"_nodb"``.

    Returns:
        tuple: A tuple of four elements.

            - allobsd (pandas.Series): ``all`` counts normalized by their own sum.
            - hitobsd (pandas.Series): ``hit`` counts normalized by ``(hit + no)``.
            - noobsd (pandas.Series): ``no`` (non-hit) counts normalized by ``(hit + no)``.
            - nodbobsd (pandas.Series): ``nodb`` counts normalized by their own sum.
    """

    header = firststep+"_gt"
    datatempgt = df[header]
    sumgt = df[header].sum()

    header = firststep+"_"+nextstep+"_all"
    datatempall = df[header]
    sumall = df[header].sum()

    header = firststep+"_"+nextstep+"_hit"
    datatemphit = df[header]
    sumhit = df[header].sum()

    header = firststep+"_"+nextstep+"_no"
    datatempno = df[header]
    sumno = df[header].sum()

    header = firststep+"_"+nextstep+"_nodb"
    datatempnodb = df[header]
    sumnodb = df[header].sum()

    gtobsd = datatempgt/sumgt
    allobsd = datatempall/sumall
    hitobsd = datatemphit/sumhit
    noobsd = datatempno/sumno
    nodbobsd = datatempnodb/sumnodb

    return(allobsd, hitobsd, noobsd, nodbobsd, gtobsd)

def generatetable(methods, totalcounter, hittotalcounter, \
                  nonhittotalcounter, nodbtotalcounter, gtcounter):
    """Flatten the per-method rank-count dictionaries into a CSV-ready table.

    Builds the header row (``"<method>_<method1>_<all|hit|no|nodb>"`` for every
    method pair) and, for ranks 1 to 50, one data row per rank containing the
    counts stored in ``totalcounter``, ``hittotalcounter``, ``nonhittotalcounter``,
    and ``nodbtotalcounter``. This is the counterpart of the columns consumed by
    :func:`obsdsetup`.

    Args:
        methods (list[str]): List of annotation-method names.
        totalcounter (dict): ``totalcounter[method][method1][rank]`` -> count of all
            observations, as produced by :func:`inputtotalcounter`.
        hittotalcounter (dict): Same structure as ``totalcounter`` but restricted to
            correct-answer (hit) observations, as produced by :func:`inputhitcounters`.
        nonhittotalcounter (dict): Same structure as ``totalcounter`` but restricted to
            incorrect-answer (non-hit) observations, as produced by
            :func:`inputhitcounters`.
        nodbtotalcounter (dict): Same structure as ``totalcounter`` but restricted to
            database-absent observations, as produced by :func:`inputtotalcounter`.

    Returns:
        list[list]: A two-dimensional list. The first row is the header; each
            subsequent row (ranks 1-50) holds the corresponding counts for every
            method pair and every category (``all``, ``hit``, ``no``, ``nodb``).
    """
    totaltext = []

    textmethod = ["rank"]
    for method in methods:
        textmethod.append(method+"_gt")

    for method in methods:
        for additonal in ["all","hit","no","nodb"]:
            for method1 in methods:
                textmethod.append(method+"_"+method1+"_"+additonal)
    totaltext = [textmethod]
    for rank1 in range(1,51):
        text = [rank1]
        for method in methods:
            text.append(gtcounter[method][rank1])

        for method in methods:
            for method1 in methods:
                text.append(totalcounter[method][method1][rank1])
            for method1 in methods:
                text.append(hittotalcounter[method][method1][rank1])
            for method1 in methods:
                text.append(nonhittotalcounter[method][method1][rank1])
            for method1 in methods:
                text.append(nodbtotalcounter[method][method1][rank1])

        #print(text)
        totaltext.append(text)
    return(totaltext)

def inputtotalcounter(method, totalcounter, thisinchikeystem, alldict):
    """Accumulate a single observation into the per-method-pair rank histogram.

    For the compound (or formula/class) that ranked 1st with ``method``, looks up
    its rank in every other method's ranking dictionary (``alldict``) and increments
    the corresponding cell of ``totalcounter``. Ranks of 50 or higher are clamped to
    50 so that the histogram has a fixed size.

    Args:
        method (str): Name of the method whose top-ranked candidate is being recorded.
        totalcounter (dict): ``totalcounter[method][targetmethod][rank]`` -> running
            count, mutated in place and also returned.
        thisinchikeystem (str): Identifier (short InChIKey, formula, or class) of the
            top-ranked candidate for ``method``.
        alldict (dict): ``alldict[targetmethod][identifier]`` -> rank, as produced by
            :func:`generatecounters`.

    Returns:
        dict: The updated ``totalcounter`` dictionary (same object as the input).
    """
    for targetmethod in alldict:
        if thisinchikeystem in alldict[targetmethod]:
            targetrank = alldict[targetmethod][thisinchikeystem]
            if targetrank >= 50:
                targetrank = 50
            totalcounter[method][targetmethod][targetrank] += 1
    return(totalcounter)


def inputhitcounters(method, hittotalcounter,nonhittotalcounter, thisinchikeystem,grandtruthinchikey_set, alldict):
    """Accumulate a single observation into the hit / non-hit rank histograms.

    Same lookup logic as :func:`inputtotalcounter`, but the count is routed to
    ``hittotalcounter`` when the top-ranked candidate matches the ground truth
    (``grandtruthinchikey_set``) and to ``nonhittotalcounter`` otherwise. This lets
    the ranking-quality statistics of a method be split into "was the annotation
    actually correct" and "was the annotation wrong" cases.

    Args:
        method (str): Name of the method whose top-ranked candidate is being recorded.
        hittotalcounter (dict): ``hittotalcounter[method][targetmethod][rank]`` ->
            running count for correct annotations, mutated in place and returned.
        nonhittotalcounter (dict): Same structure as ``hittotalcounter`` but for
            incorrect annotations, mutated in place and returned.
        thisinchikeystem (str): Identifier (short InChIKey, formula, or class) of the
            top-ranked candidate for ``method``.
        grandtruthinchikey_set (set[str]): Set of identifiers considered ground truth
            for this record (e.g. derived from a high-confidence spectral match).
        alldict (dict): ``alldict[targetmethod][identifier]`` -> rank, as produced by
            :func:`generatecounters`.

    Returns:
        tuple: A tuple of two elements.

            - hittotalcounter (dict): The updated hit-count dictionary.
            - nonhittotalcounter (dict): The updated non-hit-count dictionary.
    """
    for targetmethod in alldict:
        if thisinchikeystem in alldict[targetmethod]:
            targetrank = alldict[targetmethod][thisinchikeystem]
            if targetrank >= 50:
                targetrank = 50
            #
            # hit when top-ranked in first method is grand truth
            #
            if thisinchikeystem in grandtruthinchikey_set:
                hittotalcounter[method][targetmethod][targetrank] += 1
            #
            # non hit when top-ranked in first method is not grand truth
            #
            else:
                nonhittotalcounter[method][targetmethod][targetrank] += 1
        #
        # Here is important
        #
        #else:
        #    targetrank = 50
        #    #
        #    # hit when top-ranked in first method is grand truth
        #    #
        #    if thisinchikeystem == grandtruthinchikey:
        #        hittotalcounter[method][targetmethod][targetrank] += 1
        #    #
        #    # non hit when top-ranked in first method is not grand truth
        #    #
        #    else:
        #        nonhittotalcounter[method][targetmethod][targetrank] += 1

    return(hittotalcounter,nonhittotalcounter)

def generatecounters(methods, anno, grandtruthinchikey_set):
    """Build per-method rank lookup tables from a record's parsed annotation data.

    For each method, walks its ranked candidate list (``anno["methodrank"]``) and
    builds a "de-duplicated rank" (i.e. the rank after collapsing repeated
    identifiers, since the same identifier can legitimately appear more than once
    in the raw ranking) mapping from identifier to rank, plus the top (rank-1)
    identifier. A second, parallel set of tables (``*_nodb``) is built while
    skipping any identifier that is already known to be correct
    (``grandtruthinchikey_set``); this "no-database-hit" view is used to study
    ranking behaviour on entries that are known *not* to be the answer.

    Args:
        methods (list[str]): List of annotation-method names to process.
        anno (dict): Structured annotation data for one record, as returned by
            ``ddatoolbox.massbankio.MassBankRecord.get_proposed()``. Must contain
            the ``"methodrank"`` key: ``method -> {rank: {"tg": identifier, ...}}``.
        grandtruthinchikey_set (set[str]): Set of identifiers considered ground
            truth for this record.

    Returns:
        tuple: A tuple of four elements.

            - alldict (dict): ``alldict[method][identifier]`` -> de-duplicated rank,
              built from the full candidate list.
            - topdict (dict): ``topdict[method]`` -> identifier ranked 1st (only set
              for methods with at least one candidate).
            - alldict_nodb (dict): Same as ``alldict`` but excluding identifiers in
              ``grandtruthinchikey_set``.
            - topdict_nodb (dict): Same as ``topdict`` but computed from
              ``alldict_nodb``; defaults to ``""`` for methods with no candidates.
    """
    #
    # Counters
    #
    topdict = {}
    alldict = {}
    topdict_nodb = {}
    alldict_nodb = {}
    for method in methods:
        topdict_nodb[method] = ""
        alldict[method] = {}
        alldict_nodb[method] = {}
    #
    # Collect ranking data
    #
    for method in methods:
        counter_nodb = set()
        counter = set()
        if not method in anno["methodrank"]:
            continue
        #for rank in anno["methodrank"][method]:
        for rank in range(1,1000):
            if not rank in anno["methodrank"][method]:
                continue
            inchikey = anno["methodrank"][method][rank]["tg"]
            if re.search(r'^[A-Z]{14}\-[A-Z]{10}\-[A-Z]$', inchikey):
                thisinchikeystem = inchikey.split("-")[0]
            else:
                thisinchikeystem = inchikey
            #
            # Regular (all-candidates) record
            #
            if thisinchikeystem in counter:
                continue
            counter.add(thisinchikeystem)
            thisrank = len(counter)
            if thisrank == 1:
                topdict[method] = thisinchikeystem
            alldict[method][thisinchikeystem] = thisrank
            #
            # For the "no database (ground truth excluded)" view
            #
            if thisinchikeystem in grandtruthinchikey_set:
                continue
            if thisinchikeystem in counter_nodb:
                continue
            counter_nodb.add(thisinchikeystem)
            thisrank = len(counter_nodb)
            if thisrank == 1:
                topdict_nodb[method] = thisinchikeystem
            alldict_nodb[method][thisinchikeystem] = thisrank
    return(alldict, topdict, alldict_nodb, topdict_nodb)


def generategtpositions(methods, gtcounter, anno, grandtruthinchikey_set):
    """Build per-method rank lookup tables from a record's parsed annotation data.

    """
    #
    # Counters
    #
    for method in methods:
        counter = set()
        if not method in anno["methodrank"]:
            continue
        #for rank in anno["methodrank"][method]:
        for rank in range(1,1000):
            if not rank in anno["methodrank"][method]:
                continue
            inchikey = anno["methodrank"][method][rank]["tg"]
            if re.search(r'^[A-Z]{14}\-[A-Z]{10}\-[A-Z]$', inchikey):
                thisinchikeystem = inchikey.split("-")[0]
            else:
                thisinchikeystem = inchikey
            # Regular (all-candidates) record
            counter.add(thisinchikeystem)
            thisrank = len(counter)
            # For the "no database (ground truth excluded)" view
            if thisinchikeystem in grandtruthinchikey_set:
                if thisrank > 50:
                    thisrank = 50
                #print(method, thisrank)
                gtcounter[method][thisrank] += 1
                continue
    return(gtcounter)

def qualifybysirius(anno, filtering_method, filtering_confidence_threshold):
    """Check whether a record's top SIRIUS6 confidence score passes a threshold.

    Reads the rank-1 score of the ``"SIRIUS6_Confidence:biodatabase"`` method from
    the parsed annotation data and compares it against ``filtering_confidence_threshold``.
    This score-based filter is used to exclude low-confidence records before
    Bayesian ranking aggregation.

    Args:
        anno (dict): Structured annotation data for one record, as returned by
            ``ddatoolbox.massbankio.MassBankRecord.get_proposed()``. Must contain
            the ``"methodrank"`` key.
        filtering_method (str): Filtering-mode identifier. Currently unused by the
            active code path (see the commented-out ``"both"`` branch below) but
            kept for interface compatibility with callers.
        filtering_confidence_threshold (float): Minimum SIRIUS6 confidence score
            required to pass.

    Returns:
        bool: ``True`` if the top SIRIUS6 confidence score is greater than or equal
            to ``filtering_confidence_threshold`` (``0`` when the method is absent
            from the record), ``False`` otherwise.
    """
    #
    # Check top "SIRIUS6_Confidence:biodatabase" score
    #
    s6confidence = 0
    if "SIRIUS6_Confidence:biodatabase" in anno["methodrank"]:
        s6confidence = anno["methodrank"]["SIRIUS6_Confidence:biodatabase"][1]["score"]
    #if "SIRIUS6_Zodiac" in anno["methodrank"]:
    #    s6confidence = anno["methodrank"]["SIRIUS6_Zodiac"][1]["score"]
    #
    # qualify by the sirius confidence score
    #
    if s6confidence >= filtering_confidence_threshold:
        passing = True
    else:
        passing = False

    """

    if filtering_method == "both":
        if not (s6confidence >= filtering_confidence_threshold and s5confidence >= filtering_confidence_threshold):
            passing = False
    else:
        if not (s6confidence >= filtering_confidence_threshold or s5confidence >= filtering_confidence_threshold):
            passing = False
    """
    return(passing)



def qualifybyranking(methods, topdict, ranking_threshold, alldict):
    """Check whether any top-ranked candidate satisfies all per-method rank thresholds.

    Collects the set of candidate identifiers that were ranked 1st by at least one
    method (``topdict``), then, for each such candidate, verifies whether its rank
    in every method listed in ``ranking_threshold`` is at or better than (i.e. does
    not exceed) that method's allowed threshold. The record qualifies as soon as one
    candidate passes all thresholds simultaneously.

    Args:
        methods (list[str]): List of annotation-method names to draw candidates from.
        topdict (dict): ``topdict[method]`` -> identifier ranked 1st by ``method``,
            as produced by :func:`generatecounters`.
        ranking_threshold (dict): ``ranking_threshold[method]`` -> maximum allowed
            rank (inclusive) for ``method``.
        alldict (dict): ``alldict[method][identifier]`` -> rank, as produced by
            :func:`generatecounters`. A candidate missing from a required method's
            dictionary automatically fails the check for that method.

    Returns:
        bool: ``True`` if at least one candidate satisfies all thresholds in
            ``ranking_threshold``, ``False`` otherwise.
    """
    #
    # qualify by the ranking threshold
    #
    candidateinchikeys = set()
    for method in methods:
        thisinchikeystem = topdict[method]
        candidateinchikeys.add(thisinchikeystem)


    for candidateinchikey in candidateinchikeys:
        rankpasscheck = True
        for checkmethod in ranking_threshold:
            thresrank = ranking_threshold[checkmethod]
            #print(alldict[checkmethod], candidateinchikey)
            if candidateinchikey in alldict[checkmethod]:
                topprobrank = alldict[checkmethod][candidateinchikey]
                if topprobrank > thresrank:
                    rankpasscheck = False
            else:
                rankpasscheck = False

            #
        # If at least one candidateinchikey pass the ranking thresholds
        #
        if rankpasscheck == True:
            return(True)
    return(False)


def rankingcount(target, datasets, confidence_threshold, methods, mode = "each", inputstep = 3, rankinglimit = 5, reset = True):
    """Aggregate ranking statistics for multiple annotation methods from MassBank records.

    Reads MassBank record files for each dataset, builds the ground-truth set from
    high-confidence GNPS spectral-entropy matches, applies a SIRIUS-confidence (for
    ``target in {"compound", "formula"}``) or ranking-threshold (for
    ``target == "class"``) qualification filter, and aggregates rank-wise counts of
    all observations, hits, non-hits, and database-absent entries for every pair of
    annotation methods. This is the data source consumed by :func:`main` /
    :func:`obsdsetup` for the Bayesian Zipf-distribution fit.

    Args:
        target (str): Type of annotation to aggregate. One of ``"compound"``,
            ``"formula"``, or ``"class"``. Determines which ``get_proposed()`` view
            is read and which qualification logic is applied.
        datasets (list[str]): List of dataset identifiers to process.
        confidence_threshold (dict): Filtering threshold settings. Must contain a
            ``"confidence"`` sub-dict with keys ``"thres"`` (numeric threshold) and
            ``"method"`` (filtering mode, e.g. ``"both"`` or ``"s6"``), and a
            ``"rank"`` sub-dict mapping method name to maximum allowed rank (used
            only when ``target == "class"``).
        methods (list[str]): List of annotation method names to aggregate.
        mode (str, optional): Operation mode. When ``"each"``, results are logged to
            the per-dataset log file; otherwise to the shared statistics output file.
            Defaults to ``"each"``.
        inputstep (int, optional): Input data step number passed through to
            :func:`pathsetup`. Defaults to ``3``.
        rankinglimit (int, optional): Currently unused inside this function; kept for
            interface compatibility with callers. Defaults to ``5``.
        reset (bool, optional): Passed through to :func:`pathsetup` to control
            whether output directories are recreated. Defaults to ``True``.

    Returns:
        list[list]: A two-dimensional list where the first row is a header and subsequent
            rows contain rank-wise counts for ranks 1 to 50.
            Returns an empty list ``[]`` if no usable data is available.
    """
    totalcounter = {}
    for method in methods:
        totalcounter[method] = {}
        for method1 in methods:
            totalcounter[method][method1] = {}
            for rank1 in range(1,51):
                totalcounter[method][method1][rank1] = 0
    hittotalcounter = copy.deepcopy(totalcounter)
    nonhittotalcounter = copy.deepcopy(totalcounter)
    nodbtotalcounter = copy.deepcopy(totalcounter)

    gtcounter = {}
    for method in methods:
        gtcounter[method] = {}
        for rank1 in range(1,51):
            gtcounter[method][rank1] = 0

    totalcounter_total = 0
    totalcounter_fullannotation = 0
    totalcounter_poorconfidence = 0
    totalcounter_used = 0
    totalcounter_numofgrandtruth = 0

    filtering_confidence_threshold = confidence_threshold["confidence"]["thres"]
    filtering_method = confidence_threshold["confidence"]["method"]
    ranking_threshold = confidence_threshold["rank"]

    if target == "class":
        InChIkeydb = ddatoolbox.resources.load_inchikey2classyfire()

    for dataset in datasets:
        paths = pathsetup(dataset, inputstep = inputstep, reset = reset)

        datasetpath = glob.glob(paths["input"])
        totalfilenumber = len(datasetpath)
        thistotalcounter = 0
        for file in datasetpath:
            if not ".txt" in file:
                continue
            mb2 = ddatoolbox.massbankio.MassBankRecord()
            mb2.read_MassBankRecord(file)
            accession = mb2.get_accession()

            anno_compound = mb2.get_proposed("compound")
            #
            # collection of grandtruth from the spectral entropy > 0.7
            #
            grandtruthinchikey_set = set()
            if "Entropy:GNPS" in anno_compound["methodtg"]:
                for thisinchikeystem in anno_compound["methodtg"]["Entropy:GNPS"]:
                    score = anno_compound["methodtg"]["Entropy:GNPS"][thisinchikeystem]["score"]
                    name = anno_compound["methodtg"]["Entropy:GNPS"][thisinchikeystem]["name"]
                    if score > 0.7:
                        if target == "formula":
                            names = name.split("_")
                            if len(names) > 1:
                                thisformula = names[-1]
                                grandtruthinchikey_set.add(thisformula)
                        elif target == "compound":
                            grandtruthinchikey_set.add(thisinchikeystem)
                        elif target == "class":
                             if thisinchikeystem in InChIkeydb:
                                thisclassyfire = InChIkeydb[thisinchikeystem]
                                grandtruthinchikey_set.add(thisclassyfire)

            if len(grandtruthinchikey_set) >= 1:
                totalcounter_numofgrandtruth += 1

            anno = mb2.get_proposed(target=target)
            gtcounter = generategtpositions(methods, gtcounter, anno, grandtruthinchikey_set)
            #
            # Check whether the record was detected by all required methods
            #
            totalcounter_total += 1
            thistotalcounter += 1
            checker = []
            for method in methods:
                if not method in anno["methodtg"]:
                    checker.append(method)
            if len(checker) > 0:
                totalcounter_fullannotation += 1
                continue
            #
            # Counters
            #
            alldict, topdict, alldict_nodb, topdict_nodb = generatecounters(methods, anno, grandtruthinchikey_set)

            #
            # Check top "SIRIUS6_Confidence:biodatabase" score
            #
            if target == "formula" or target == "compound":
                #
                # SIRIUS 6 confidence score is used for qualification of formulan and compund
                #
                if qualifybysirius(anno_compound, filtering_method, filtering_confidence_threshold) == False:
                    #print(filtering_method, filtering_confidence_threshold)
                    totalcounter_poorconfidence += 1
                    continue
            if target == "class":
                if qualifybyranking(methods, topdict, ranking_threshold, alldict) == False:
                    continue

            totalcounter_used += 1
            print("DDAZipf_{:s}: Chech ranking ({:d}/{:d})".format(target, thistotalcounter, totalfilenumber), accession, grandtruthinchikey_set)
            #
            # Check the rank for each method
            #
            for method in methods:
                thisinchikeystem = topdict[method]
                if thisinchikeystem == "":
                    continue
                # All data
                totalcounter = inputtotalcounter(method, totalcounter, thisinchikeystem, alldict)

                #for targetmethod in alldict:
                #    if thisinchikeystem in alldict[targetmethod]:
                #        targetrank = alldict[targetmethod][thisinchikeystem]
                #        if targetrank >= 50:
                #            targetrank = 50
                #        if method == targetmethod and targetrank != 1:
                #            print(accession, method,targetmethod,targetrank)

                # Only when there is a ground truth available
                if len(grandtruthinchikey_set) == 0:
                    continue
                hittotalcounter,nonhittotalcounter = \
                inputhitcounters(method, hittotalcounter,nonhittotalcounter, thisinchikeystem,grandtruthinchikey_set, alldict)
                # Top-ranked inchikey in top5dict_nodb
                thisinchikeystem = topdict_nodb[method]
                if thisinchikeystem == "":
                    continue
                nodbtotalcounter = inputtotalcounter(method, nodbtotalcounter, thisinchikeystem, alldict_nodb)

    if totalcounter_used == 0:
        return([])

    totaltext = generatetable(methods, totalcounter, hittotalcounter, nonhittotalcounter, nodbtotalcounter, gtcounter)
    
    # Writing to a CSV file
    print(dataset, confidence_threshold,  totalcounter_total, totalcounter_numofgrandtruth, totalcounter_fullannotation, totalcounter_poorconfidence, totalcounter_used)
    if mode == "each":
        logfilename = paths["logfile"]
    else:
        logfilename = paths["statisticsoutput"]

    dt_now = datetime.datetime.now()
    ddatoolbox.ddatoolbox.writeLog(logfilename, ["ZipfRanking", str(datasets), dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "confidence_threshold", confidence_threshold,  "total", totalcounter_total, "numofgrandtruth", totalcounter_numofgrandtruth, "fullannotation", totalcounter_fullannotation, "poorconfidence",totalcounter_poorconfidence, "Used", totalcounter_used], mode = "a")

    return(totaltext)



def distoutput(paths, trace, abbreviations, df, rankinglimit):
    """Plot posterior distributions for every method pair and summarize them as a table.

    For every ordered ``(firststep, nextstep)`` method combination, draws a 2x5 grid
    of subplots per figure: log-log scatter plots of the observed ``all``/``hit``/
    ``no``/``nodb`` frequency distributions (see :func:`obsdsetup`) side by side with
    histograms of the posterior samples for the Zipf shape parameters ``la``, ``lb``,
    ``l0``, ``l1``, and ``alpha`` (see :func:`main`). Each figure is saved as a PNG
    file under ``paths["softwareoutputpath"]``. In addition, the mean and 95%
    credible interval of every sampled parameter are collected into a summary table.

    Args:
        paths (dict): Dictionary containing file path information, as produced by
            :func:`pathsetup`. Must include key ``"softwareoutputpath"``
            (``pathlib.Path``) pointing to the output directory for the PNG files.
        trace (dict): Dictionary mapping parameter names to arrays of posterior
            samples, i.e. the ``tracedict`` returned by :func:`main`.
        abbreviations (dict): Dictionary mapping method names to their abbreviations.
        df (pandas.DataFrame): DataFrame containing ranking-count data, in the same
            format expected by :func:`obsdsetup`.
        rankinglimit (int): Upper rank limit (number of leading points) used for the
            scatter plots of observed data.

    Returns:
        list[list]: A two-dimensional list where each row has the format
            ``[firststep, nextstep, param_name, mean, 2.5th_percentile, 97.5th_percentile]``
            for each of the five parameters (``la``, ``lb``, ``l0``, ``l1``, ``alpha``)
            of every method combination.
    """
    outputtext = []
    allobsd = {}
    hitobsd = {}
    noobsd = {}
    nodbobsd = {}
    gtobsd = {}
    methods = list(abbreviations.keys())
    combinations = list(itertools.permutations(methods, 2))
    for firststep, nextstep in combinations:
        print(firststep, nextstep)
        cmb = (firststep, nextstep)
        allobsd[cmb], hitobsd[cmb], noobsd[cmb], nodbobsd[cmb], gtobsd[firststep] = obsdsetup(firststep, nextstep, df)

    for firststep, nextstep in combinations:
        cmb = (firststep, nextstep)
        firststepabr = abbreviations[firststep]
        cmbabr = abbreviations[firststep]+abbreviations[nextstep]

        ax = plt.subplot(251)
        xaxis =np.log(range(1,rankinglimit+1))
        plt.scatter(xaxis, np.log(allobsd[cmb][:rankinglimit]))

        plt.title("all")

        ax = plt.subplot(252)
        plt.scatter(xaxis, np.log(hitobsd[cmb][:rankinglimit]))
        #print(cmbabr)
        #print(xaxis)
        #print(hitobsd[cmb][:rankinglimit])
        plt.title("hit")

        ax = plt.subplot(253)
        plt.scatter(xaxis, np.log(noobsd[cmb][:rankinglimit]))
        plt.title("no")

        ax = plt.subplot(254)
        plt.scatter(xaxis, np.log(nodbobsd[cmb][:rankinglimit]))
        plt.title("nodb")

        ax = plt.subplot(255)
        targetparam = 'la_'+firststepabr
        alpha_samples = trace[targetparam]

        tempdata = [firststep, nextstep, targetparam,np.mean(alpha_samples), np.percentile(alpha_samples, 2.5), np.percentile(alpha_samples, 97.5)]
        print(tempdata)
        outputtext.append(tempdata)

        plt.hist(alpha_samples, histtype='stepfilled', bins=30, alpha=0.85, label=targetparam, color="#A60628", density=True)
        plt.title(targetparam)

        ax = plt.subplot(256)
        targetparam = 'lb_'+cmbabr
        alpha_samples = trace[targetparam]
        tempdata = [firststep, nextstep, targetparam,np.mean(alpha_samples), np.percentile(alpha_samples, 2.5), np.percentile(alpha_samples, 97.5)]
        print(tempdata)
        outputtext.append(tempdata)

        plt.hist(alpha_samples, histtype='stepfilled', bins=30, alpha=0.85, label=targetparam, color="#A60628", density=True)
        plt.title(targetparam)

        ax = plt.subplot(257)
        targetparam = 'l0_'+cmbabr
        alpha_samples = trace[targetparam]
        tempdata = [firststep, nextstep, targetparam,np.mean(alpha_samples), np.percentile(alpha_samples, 2.5), np.percentile(alpha_samples, 97.5)]
        print(tempdata)
        outputtext.append(tempdata)
        plt.hist(alpha_samples, histtype='stepfilled', bins=30, alpha=0.85, label=targetparam, color="#A60628", density=True)
        plt.title(targetparam)

        ax = plt.subplot(258)
        targetparam = 'l1_'+cmbabr
        alpha_samples = trace[targetparam]
        tempdata = [firststep, nextstep, targetparam,np.mean(alpha_samples), np.percentile(alpha_samples, 2.5), np.percentile(alpha_samples, 97.5)]
        print(tempdata)
        outputtext.append(tempdata)
        plt.hist(alpha_samples, histtype='stepfilled', bins=30, alpha=0.85, label=targetparam, color="#A60628", density=True)
        plt.title(targetparam)

        ax = plt.subplot(259)
        targetparam = 'alpha'
        alpha_samples = trace[targetparam]
        tempdata = [firststep, nextstep, targetparam,np.mean(alpha_samples), np.percentile(alpha_samples, 2.5), np.percentile(alpha_samples, 97.5)]
        print(tempdata)
        outputtext.append(tempdata)
        plt.hist(alpha_samples, histtype='stepfilled', bins=30, alpha=0.85, label=targetparam, color="#A60628", density=True)
        plt.title(targetparam)

        filename = cmbabr+".png"

        figoutputpath = paths["softwareoutputpath"] / filename
        plt.savefig(figoutputpath)
        plt.clf()
    return(outputtext)


def zipfsetup(datasets, abbreviations, thresholds, target = "compound", inputstep = 3, rankinglimit = 10, reset = True):
    """Generate and save Zipf posterior distributions for compound annotation in batch mode.

    For each instrument type (TOF and FT), aggregates ranking statistics from the matching
    datasets and estimates the posterior distribution via MCMC sampling with PyMC.
    Results are saved as a pickle file (tracedict) and a CSV file (posterior).

    Args:
        datasets (list[str]): List of dataset identifiers to process.
        abbreviations (dict): Dictionary mapping method names to their abbreviations.
        confidence_threshold (dict): Per-instrument-type filtering threshold dictionary.
            Keys ``"TOF"`` and ``"FT"`` each map to ``{"thres": float, "method": str}``.
        inputstep (int, optional): Input data step number. Defaults to ``3``.
        rankinglimit (int, optional): Upper rank limit used for Zipf model fitting.
            Defaults to ``10``.
        reset (bool, optional): When ``True``, skips ranking aggregation and reuses the
            existing CSV file, running only MCMC sampling. Defaults to ``True``.

    Returns:
        None
    """

    methods_compound = list(abbreviations["compound"].keys())
    methods_formula = list(abbreviations["formula"].keys())
    methods_class = list(abbreviations["class"].keys())

    threshold_compound = thresholds["compound"]
    threshold_formula = thresholds["formula"]
    threshold_class = thresholds["class"]



    respath = ddatoolbox.paths.generate_paths_resource()
    for mode in ["FT", "TOF"]:
        if mode == "TOF":
            filenameranking = respath["zipfroot"] +"/zipftof_{}/ranking.csv".format(target)
            filenametrace = respath["zipfroot"] +"/zipftof_{}/sample.pkl".format(target)
            filenameposterior = respath["zipfroot"] +"/zipftof_{}/posterior.csv".format(target)
            filenamelog = respath["zipfroot"] +"/zipftof_{}/log.txt".format(target)
            paths = {}
            paths["softwareoutputpath"] = Path(respath["zipfroot"] +"/zipftof_{}/".format(target)).resolve()
        else:
            filenameranking = respath["zipfroot"] +"/zipfft_{}/ranking.csv".format(target)
            filenametrace =respath["zipfroot"] +"/zipfft_{}/sample.pkl".format(target)
            filenameposterior =respath["zipfroot"] +"/zipfft_{}/posterior.csv".format(target)
            filenamelog = respath["zipfroot"] +"/zipfft_{}/log.txt".format(target)
            paths = {}
            paths["softwareoutputpath"] = Path(respath["zipfroot"] +"/zipfft_{}/".format(target)).resolve()

        datasets_new = []
        for dataset in datasets:
            template = ddatoolbox.massbankio.read_RecordTemplate_dataset(dataset)
            thispol = template["INSTRUMENT_TYPE"]
            if mode in thispol:
                datasets_new.append(dataset)
        if len(datasets_new) == 0:
            print("Skipped", mode)
            continue

        if reset == True:
            if target == "compound":
                totaltext = rankingcount(target, datasets_new, threshold_compound[mode], methods_compound, mode = "batch", inputstep = inputstep, rankinglimit = rankinglimit, reset = reset)
            elif target == "formula" or target == "mf":
                target = "formula"
                #totaltext = rankingcount_formula(datasets_new, threshold_formula[mode], methods_compound, methods_formula, mode = "batch", inputstep = inputstep, rankinglimit = rankinglimit, reset = reset)
                totaltext = rankingcount(target, datasets_new, threshold_formula[mode], methods_formula, mode = "batch", inputstep = inputstep, rankinglimit = rankinglimit, reset = reset)
            elif target == "class":
                #totaltext = rankingcount_class(datasets_new, threshold_class[mode], methods_compound, methods_class, mode = "batch", inputstep = inputstep, rankinglimit = rankinglimit, reset = reset)
                totaltext = rankingcount(target, datasets_new, threshold_class[mode], methods_class, mode = "batch", inputstep = inputstep, rankinglimit = rankinglimit, reset = reset)


            with open(filenameranking, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerows(totaltext)

        df = pd.read_csv(filenameranking)

        trace, tracedict = main(abbreviations[target], df, rankinglimit, target = target)
        with open(filenametrace, 'wb') as f:
            pickle.dump(tracedict, f)

        with open(filenametrace, 'rb') as f:
            tracedict = pickle.load(f)
        posterioroutput_text = distoutput(paths, tracedict, abbreviations[target], df, rankinglimit)
        with open(filenameposterior, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerows(posterioroutput_text)
        dt_now = datetime.datetime.now()
        ddatoolbox.ddatoolbox.writeLog(filenamelog, ["Created", str(datasets_new), dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "Zipf"], mode = "a")


if __name__ == '__main__':

    description = """
    ShinMassBank pipeline Zipf (generation of the preestimated Zipf distribution)\n
    Datasets must be annotated by SIRIUS 5 SIRIUS6, metfrag, MS-FINDER, and retreived. \n
    Please use the preestimated Zipf distribution for a normal use by -a option of DDAcalcfdr.py.
    Example 2: Generation of a new preestimated Zipf distribution \n
    python DDAZipf.py HSA001 HSA002 -p\n
    \n
    """
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("datasets", nargs='*', default=[])    # positional argument: dataset identifiers
    parser.add_argument("-p", "--preestimated", action='store_true', help='Generation of a new preestimated Zipf distribution')    # required option

    args = parser.parse_args()    # Parse command-line arguments

    datasets = args.datasets

    parameter = setparameter()

    thresholds = {}
    thresholds["compound"] = parameter["confidence_threshold"]
    thresholds["formula"] = parameter["formulaconfidence_threshold"]
    thresholds["class"] = parameter["class_threshold"]

    rankinglimit = parameter["rankinglimit"]
    rankinglimit_formula = parameter["rankinglimit_formula"]
    rankinglimit_class = parameter["rankinglimit_class"]

    inputstep = parameter["inputstep"]

    abbreviations = {}
    abbreviations["compound"] = parameter["abbreviations"]
    abbreviations["formula"] = parameter["abbreviations_formula"]
    abbreviations["class"] = parameter["abbreviations_class"]

    if args.preestimated == True:

        zipfsetup(datasets, abbreviations, thresholds, target = "class", inputstep = inputstep, rankinglimit = rankinglimit_class, reset = False)
        zipfsetup(datasets, abbreviations, thresholds, target = "compound", inputstep = inputstep, rankinglimit = rankinglimit, reset = False)
        zipfsetup(datasets, abbreviations, thresholds, target = "formula", inputstep = inputstep, rankinglimit = rankinglimit_formula, reset = False)
