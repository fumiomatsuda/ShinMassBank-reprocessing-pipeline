
#-------------------------------------------------------------------------------
# Name:        compoundsearch.py
# Purpose:
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import glob, csv, pickle, time, sys
import subprocess
from pathlib import Path
import os, shutil, re, datetime, math
from scipy.stats import spearmanr, norm
import joblib

from . import massbankio
from . import mfsearch
from . import ddatoolbox
from . import massbankserver
from . import paths

rankingthreshold = 100
production_correlation_alpha = 0.01

def login_sirius(version, address, path):
    """Log into the SIRIUS software (version 5 or 6) using the provided credentials.

    Constructs and executes the SIRIUS login command via the command line,
    using the specified user email address and the appropriate executable
    path for the given SIRIUS version.

    Args:
        version (int): The version of SIRIUS to use (5 or 6).
        address (str): User email address for login.
        path (dict): Dictionary containing paths to SIRIUS executables.
            Must include keys ``"sirius5"`` and ``"sirius6"`` as pathlib.Path objects.

    Returns:
        None
    """
    #'"C:\Program Files\sirius\sirius" login -u fmatsuda@ist.osaka-u.ac.jp -p'
    if version == 6:
        exec = path["sirius6"]/"sirius"
    else:
        exec = path["sirius5"]/"sirius"

    commanttext = '"{exec:s}" login -u {address:s} -p'
    commanttext = commanttext.format(exec = str(exec), address = address)

    print(commanttext)
    subprocess.run(commanttext, shell=True)



def exec_sirius(version, dataset, path):
    """Configure and execute SIRIUS (version 5 or 6) for formula and structure annotation.

    Reads instrument type, polarity, and m/z thresholds from the RecordTemplate and
    input MassBank records, constructs the appropriate SIRIUS configuration command
    for the given version and ionization mode, and runs both the configuration and
    analysis pipeline via subprocess. Zodiac is included when the dataset contains
    100 or more records.

    Args:
        version (int): SIRIUS version to use (5 or 6).
        dataset (str): Dataset identifier string.
        path (dict): Dictionary containing file path information with the following keys:
            ``"dataset"`` (pathlib.Path), ``"software"`` (pathlib.Path),
            ``"logfilename"`` (str), ``"startlogfilename"`` (str),
            ``"output"`` (pathlib.Path), ``"softwareinput"`` (pathlib.Path),
            ``"softwareoutput"`` (pathlib.Path), ``"sirius5"`` (pathlib.Path),
            ``"sirius6"`` (pathlib.Path).

    Returns:
        None
    """
    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    softwareinput = path["softwareinput"]
    softwareoutput = path["softwareoutput"]

    softwareoutputfile = path["softwareoutput"]/"output"

    if version == 6:
        exec = path["sirius6"]/"sirius"
    else:
        exec = path["sirius5"]/"sirius"

    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")
    #
    # Check polarity
    #

    files = glob.glob(str(inputpath))
    machine = "orbitrap"



    for file in files:
        if not ".txt" in file:
            continue
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(file)
        #
        # Getting information
        #
        polarity = mb2.get_polarity()
        #ms2thres = float(mb2.get_anything("COMMENT:","MS2_M/Z_AVERAGED_THRES"))
        #ms1thres = float(mb2.get_anything("COMMENT:","MS1_M/Z_AVERAGED_THRES"))
        if "QTOF" in mb2.get_instrument_type():
            machine = "qtof"

        break
    ms1thres = float(settings["ms1_mz_averaged_thres"])
    ms2thres = float(settings["ms2_mz_averaged_thres"])
    #
    # If dataset has more than 100 records, Zodiac is employed
    #
    numoffimes = len(files)
    zodiacexec = "zodiac "

    if numoffimes < 100:
        zodiacexec = ""


    if version == 6:
        #config --FormulaSearchSettings.applyFormulaConstraintsToBottomUp=false --IsotopeSettings.filter=false --UseHeuristic.useOnlyHeuristicAboveMz=650 --FormulaSearchDB=, --Timeout.secondsPerTree=0 --FormulaSettings.enforced=H,C,N,O,P --Timeout.secondsPerInstance=0 --AlgorithmProfile=orbitrap --SpectralMatchingMassDeviation.allowedPeakDeviation=5.0ppm --AdductSettings.ignoreDetectedAdducts=false --AdductSettings.enforced=, --AdductSettings.prioritizeInputFileAdducts=true --UseHeuristic.useHeuristicAboveMz=300 --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation=5.0ppm --SpectralMatchingMassDeviation.allowedPrecursorDeviation=5.0ppm --FormulaSearchSettings.performDeNovoBelowMz=400.0 --FormulaSearchSettings.applyFormulaConstraintsToDatabaseCandidates=false --EnforceElGordoFormula=true --NumberOfCandidatesPerIonization=1 --FormulaSettings.detectable=B,S,Cl,Se,Br --NumberOfCandidates=10 --AdductSettings.fallback=[[M+H]+] --FormulaSearchSettings.performBottomUpAboveMz=0 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --FormulaResultThreshold=true --ExpansiveSearchConfidenceMode.confidenceScoreSimilarityMode=APPROXIMATE --StructureSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB --RecomputeResults=false
        #siriusini_pos = '"{exec:s}"  config --FormulaSearchSettings.applyFormulaConstraintsToBottomUp=false --IsotopeSettings.filter=true --UseHeuristic.useOnlyHeuristicAboveMz=650 --AdductSettings.detectable=[[M+H]+,[M-H2O+H]+,[M-H4O2+H]+,[M+H3N+H]+] --FormulaSearchDB=, --Timeout.secondsPerTree=0 --FormulaSettings.enforced=H,C,N,O,P --Timeout.secondsPerInstance=0 --AlgorithmProfile={machine:s} --SpectralMatchingMassDeviation.allowedPeakDeviation={ms2thres:.2f}ppm --AdductSettings.ignoreDetectedAdducts=false --AdductSettings.enforced=, --AdductSettings.prioritizeInputFileAdducts=true --UseHeuristic.useHeuristicAboveMz=300 --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation={ms2thres:.2f}ppm --SpectralMatchingMassDeviation.allowedPrecursorDeviation={ms1thres:.2f}ppm --FormulaSearchSettings.performDeNovoBelowMz=400.0 --FormulaSearchSettings.applyFormulaConstraintsToDatabaseCandidates=false --EnforceElGordoFormula=true --NumberOfCandidatesPerIonization=1 --FormulaSettings.detectable=B,S,Cl,Se,Br --NumberOfCandidates=10 --AdductSettings.fallback=[[M+H]+] --FormulaSearchSettings.performBottomUpAboveMz=0 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --FormulaResultThreshold=true --ExpansiveSearchConfidenceMode.confidenceScoreSimilarityMode=APPROXIMATE --StructureSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB --RecomputeResults=true'
        #siriusini_neg = '"{exec:s}"  config --FormulaSearchSettings.applyFormulaConstraintsToBottomUp=false --IsotopeSettings.filter=true --UseHeuristic.useOnlyHeuristicAboveMz=650 ---AdductSettings.detectable=[[M-H]-,[M-H2O-H]-] -FormulaSearchDB=, --Timeout.secondsPerTree=0 --FormulaSettings.enforced=H,C,N,O,P --Timeout.secondsPerInstance=0 --AlgorithmProfile={machine:s} --SpectralMatchingMassDeviation.allowedPeakDeviation={ms2thres:.2f}ppm --AdductSettings.ignoreDetectedAdducts=false --AdductSettings.enforced=, --AdductSettings.prioritizeInputFileAdducts=true --UseHeuristic.useHeuristicAboveMz=300 --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation={ms2thres:.2f}ppm --SpectralMatchingMassDeviation.allowedPrecursorDeviation={ms1thres:.2f}ppm --FormulaSearchSettings.performDeNovoBelowMz=400.0 --FormulaSearchSettings.applyFormulaConstraintsToDatabaseCandidates=false --EnforceElGordoFormula=true --NumberOfCandidatesPerIonization=1 --FormulaSettings.detectable=B,S,Cl,Se,Br --NumberOfCandidates=10 --AdductSettings.fallback=[[M-H]-] --FormulaSearchSettings.performBottomUpAboveMz=0 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --FormulaResultThreshold=true --ExpansiveSearchConfidenceMode.confidenceScoreSimilarityMode=APPROXIMATE --StructureSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB --RecomputeResults=true'
        #spectra-search formulas zodiac fingerprints classes structures
        #exectext = '"{exec:s}" --recompute --input {input:s} --output {outputfile:s} formula {zodiac:s} fingerprint classes structures summaries -o {output:s} --top-k-summary=100'

        #config --AlgorithmProfile=qtof --MS2MassDeviation.allowedMassDeviation=10.0ppm --SpectralSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB --AdductSettings.fallback=[[M+H]+,[M+Na]+,[M+K]+] --FormulaSettings.enforced=H,C,N,O,P --IdentitySearchSettings.precursorDeviation=20.0ppm --FormulaSearchSettings.performBottomUpAboveMz=0 --FormulaSearchDB=, --StructureSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB formulas fingerprints classes structures
        #siriusini_pos = '"{exec:s}"  config --AlgorithmProfile={machine:s} --MS2MassDeviation.allowedMassDeviation={ms2thres:.2f}ppm --SpectralSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB --AdductSettings.fallback=[[M+H]+,[M+Na]+,[M+K]+] --FormulaSettings.enforced=H,C,N,O,P --RecomputeResults=true --IdentitySearchSettings.precursorDeviation={ms1thres:.2f}ppm --FormulaSearchSettings.performBottomUpAboveMz=0 --FormulaSearchDB=, --StructureSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB'

        #config --AlgorithmProfile=qtof --MS2MassDeviation.allowedMassDeviation=10.0ppm --SpectralSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB --AdductSettings.fallback=[[M-H]-] --AdductSettings.enforced=[[M-H]-] --FormulaSettings.enforced=H,C,N,O,P --RecomputeResults=true --IdentitySearchSettings.precursorDeviation=20.0ppm --FormulaSearchSettings.performBottomUpAboveMz=0 --FormulaSearchDB=, --StructureSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB formulas zodiac fingerprints classes structures
        #siriusini_neg = '"{exec:s}"  config --AlgorithmProfile={machine:s} --MS2MassDeviation.allowedMassDeviation={ms2thres:.2f}ppm --SpectralSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB --AdductSettings.fallback=[[M-H]-] --AdductSettings.enforced=[[M-H]-] --FormulaSettings.enforced=H,C,N,O,P --RecomputeResults=true --IdentitySearchSettings.precursorDeviation={ms1thres:.2f}ppm --FormulaSearchSettings.performBottomUpAboveMz=0 --FormulaSearchDB=, --StructureSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,SUPERNATURAL,TeroMol,YMDB'


        exectext = paths.getsettings("RESOURCES", "exectext")
        siriusini_pos = paths.getsettings("RESOURCES", "siriusini_pos")
        siriusini_neg = paths.getsettings("RESOURCES", "siriusini_neg")
        #exectext= str(Path(exectext).resolve())
        #siriusini_pos= str(Path(siriusini_pos).resolve())
        #siriusini_neg= str(Path(siriusini_neg).resolve())

    else:

        #siriusini_pos = '"C:\Program Files\sirius\sirius"  config --IsotopeSettings.filter=true --FormulaSearchDB=BIO --Timeout.secondsPerTree=0 --FormulaSettings.enforced=HCNOP --Timeout.secondsPerInstance=0 --AdductSettings.detectable=[[M+H]+,[M-H2O+H]+,[M-H4O2+H]+,[M+H3N+H]+] --UseHeuristic.mzToUseHeuristicOnly=650 --AlgorithmProfile=qtof --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation=10.0ppm --NumberOfCandidatesPerIon=1 --UseHeuristic.mzToUseHeuristic=300 --FormulaSettings.detectable=B,Cl,Br,Se,S --NumberOfCandidates=20 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacEpochs.burnInPeriod=2000 --ZodiacEpochs.numberOfMarkovChains=10 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --ZodiacEpochs.iterations=20000 --AdductSettings.enforced=, --AdductSettings.fallback=[[M+H]+] --FormulaResultThreshold=true --InjectElGordoCompounds=true --StructureSearchDB=BIO --RecomputeResults=True'
        #siriusini_neg = '"C:\Program Files\sirius\sirius"  config --IsotopeSettings.filter=true --FormulaSearchDB=BIO --Timeout.secondsPerTree=0 --FormulaSettings.enforced=HCNOP --Timeout.secondsPerInstance=0 --AdductSettings.detectable=[[M-H]-,[M-H2O-H]-] --UseHeuristic.mzToUseHeuristicOnly=650 --AlgorithmProfile=qtof --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation=10.0ppm --NumberOfCandidatesPerIon=1 --UseHeuristic.mzToUseHeuristic=300 --FormulaSettings.detectable=B,Cl,Br,Se,S --NumberOfCandidates=20 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacEpochs.burnInPeriod=2000 --ZodiacEpochs.numberOfMarkovChains=10 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --ZodiacEpochs.iterations=20000 --AdductSettings.enforced=, --AdductSettings.fallback=[[M-H]-] --FormulaResultThreshold=true --InjectElGordoCompounds=true --StructureSearchDB=BIO --RecomputeResults=True'

        siriusini_pos = '"{exec:s}"  config --IsotopeSettings.filter=true --FormulaSearchDB=BIO --Timeout.secondsPerTree=0 --FormulaSettings.enforced=HCNOP --Timeout.secondsPerInstance=0 --AdductSettings.detectable=[[M+H]+,[M-H2O+H]+,[M-H4O2+H]+,[M+H3N+H]+] --UseHeuristic.mzToUseHeuristicOnly=650 --AlgorithmProfile={machine:s} --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation={ms2thres:.2f}ppm --NumberOfCandidatesPerIon=1 --UseHeuristic.mzToUseHeuristic=300 --FormulaSettings.detectable=B,Cl,Br,Se,S --NumberOfCandidates=20 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacEpochs.burnInPeriod=2000 --ZodiacEpochs.numberOfMarkovChains=10 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --ZodiacEpochs.iterations=20000 --AdductSettings.enforced=, --AdductSettings.fallback=[[M+H]+] --FormulaResultThreshold=true --InjectElGordoCompounds=true --StructureSearchDB=BIO --RecomputeResults=True'
        siriusini_neg = '"{exec:s}"  config --IsotopeSettings.filter=true --FormulaSearchDB=BIO --Timeout.secondsPerTree=0 --FormulaSettings.enforced=HCNOP --Timeout.secondsPerInstance=0 --AdductSettings.detectable=[[M-H]-,[M-H2O-H]-] --UseHeuristic.mzToUseHeuristicOnly=650 --AlgorithmProfile={machine:s} --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation={ms2thres:.2f}pppm --NumberOfCandidatesPerIon=1 --UseHeuristic.mzToUseHeuristic=300 --FormulaSettings.detectable=B,Cl,Br,Se,S --NumberOfCandidates=20 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacEpochs.burnInPeriod=2000 --ZodiacEpochs.numberOfMarkovChains=10 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --ZodiacEpochs.iterations=20000 --AdductSettings.enforced=, --AdductSettings.fallback=[[M-H]-] --FormulaResultThreshold=true --InjectElGordoCompounds=true --StructureSearchDB=BIO --RecomputeResults=True'


        #                                                   config --IsotopeSettings.filter=true --FormulaSearchDB=BIO --Timeout.secondsPerTree=0 --FormulaSettings.enforced=HCNOP --Timeout.secondsPerInstance=0 --AdductSettings.detectable=[[M+H]+,[M+H3N+H]+,[M-H2O+H]+,[M-H4O2+H]+] --UseHeuristic.mzToUseHeuristicOnly=650 --AlgorithmProfile=qtof --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation=10.0ppm --NumberOfCandidatesPerIon=1 --UseHeuristic.mzToUseHeuristic=300 --FormulaSettings.detectable=B,Cl,Br,Se,S --NumberOfCandidates=20 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacEpochs.burnInPeriod=2000 --ZodiacEpochs.numberOfMarkovChains=10 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --ZodiacEpochs.iterations=20000 --AdductSettings.enforced=, --AdductSettings.fallback=[[M+H]+] --FormulaResultThreshold=true --InjectElGordoCompounds=true --StructureSearchDB=BIO --RecomputeResults=false formula zodiac fingerprint structure canopus
        #exectext_orignal = '"C:\Program Files\sirius\sirius" --recompute --input D:\MassBank\[dataset]\SIRIUS5ms\ --output D:\MassBank\[dataset]\SIRIUS5 formula  --profile qtof --database BIO --candidates 10 zodiac fingerprint structure --database BIO canopus write-summaries -o D:\MassBank\[dataset]\SIRIUS5'
        exectext = '"{exec:s}" --recompute --input {input:s} --output {output:s} formula --profile {machine:machine} --database BIO --candidates 10 {zodiac:s} fingerprint structure --database BIO canopus write-summaries -o {output:s}'


    #
    # SetUp SIRIUS for target
    #

    if polarity == "POSITIVE":
        siriusini = str(siriusini_pos)

    else:
        siriusini = str(siriusini_neg)



    siriusini = siriusini.format(exec = str(exec), ms2thres = ms2thres, ms1thres = ms1thres, machine = machine)
    exectext = exectext.format(zodiac = str(zodiacexec), machine = machine, exec = str(exec), input = str(softwareinput), output = str(softwareoutput), outputfile = str(softwareoutputfile))


    print(siriusini)
    print(exectext)

    subprocess.run(siriusini, shell=True)
    subprocess.run(exectext, shell=True)




def output_sirius6(dataset, path):
    """Parse SIRIUS 6 output files and write compound/formula annotations to output.txt.

    Reads the top-100 structure and formula identification TSV files produced by SIRIUS 6,
    extracts Confidence scores, CSI:FingerID scores, Zodiac scores, and Sirius scores,
    and appends all ranked results as tab-separated rows to the dataset's ``output.txt`` file.

    Args:
        dataset (str): Dataset identifier string.
        path (dict): Dictionary containing file path information with the following keys:
            ``"dataset"`` (pathlib.Path), ``"software"`` (pathlib.Path),
            ``"logfilename"`` (str), ``"startlogfilename"`` (str),
            ``"output"`` (pathlib.Path), ``"softwareinput"`` (pathlib.Path),
            ``"softwareoutput"`` (pathlib.Path).

    Returns:
        None
    """
    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    softwareinput = path["softwareinput"]
    softwareoutput = path["softwareoutput"]
    #
    #
    # SetUp MSFINDER.INI
    #


    accession_stem_ms2 = dataset+"P"
    accession_stem_ms1 = dataset+"S"


    confidencedict = {}
    csidict = {}

    zodiacdict = {}
    siriusdict = {}

    outputtable = []



    #print(paths, softwareoutput)





    structure_candidates_path = softwareoutput / "structure_identifications_top-100.tsv"
    #
    # If there is no result
    #
    if not os.path.isfile(structure_candidates_path):
        print("No output", structure_candidates_path)
        return()
    #
    # Read structure_candidates.tsv
    #
    rows = ddatoolbox.readtable(structure_candidates_path,\
     ["structurePerIdRank","molecularFormula", "InChIkey2D","ConfidenceScoreApproximate", \
     "adduct","CSI:FingerIDScore","name","mappingFeatureId"], filetype = "tab")


    for row in rows:
        if int(row["structurePerIdRank"]) == 1:
            confidence=float(row["ConfidenceScoreApproximate"])
        else:
            confidence = 0.0

        rank = int(row["structurePerIdRank"])
        csi = float(row["CSI:FingerIDScore"])
        inchikey = row["InChIkey2D"]
        Inchikey1_query = inchikey.split("-")[0]

        ionform=row["adduct"].replace(" ", "")
        formula=row["molecularFormula"]
        name = row["name"]
        mappingFeatureId = row["mappingFeatureId"]
        # 655215408016792237_MSBNK-SMB_Measured-HSA001P0205016_1_MSBNK-SMB_Measured-HSA001P0205016_1
        accession = re.findall(r'(MSBNK-[A-Za-z0-9_]+-[A-Za-z0-9]+)', mappingFeatureId)
        accession = accession[0]
        #print(accession, confidence, row["ConfidenceScoreApproximate"], row["structurePerIdRank"])

        if confidence > 0.0:
            if not accession in confidencedict:
                confidencedict[accession] = []
            temp = {"score":confidence, "ionform": ionform, "name": name, "formula": formula, "InChIKey": Inchikey1_query}
            confidencedict[accession].append(temp)

        if rank < rankingthreshold:
            if not accession in csidict:
                csidict[accession] = []
            temp = {"score":csi, "ionform": ionform, "name": name, "formula": formula, "InChIKey": Inchikey1_query}
            csidict[accession].append(temp)

    #
    # mf
    #
    formula_candidates_path = softwareoutput / "formula_identifications_top-100.tsv"
    #
    # If there is no result
    #
    if not os.path.isfile(formula_candidates_path):
        print("No output", formula_candidates_path)
        return()
    #
    # Read formula_identifications_top-100.tsv
    #
    rows = ddatoolbox.readtable(formula_candidates_path,\
     ["formulaRank","molecularFormula", "adduct", "SiriusScore","ZodiacScore","mappingFeatureId"])


    for row in rows:

        if row["ZodiacScore"] == "nan":
            zodiac=0.0
        else:
            zodiac=float(row["ZodiacScore"])
        rank = int(row["formulaRank"])
        sirius = float(row["SiriusScore"])
        ionform=row["adduct"].replace(" ", "")
        formula=row["molecularFormula"]
        mappingFeatureId = row["mappingFeatureId"]
        # 655215408016792237_MSBNK-SMB_Measured-HSA001P0205016_1_MSBNK-SMB_Measured-HSA001P0205016_1
        accession = re.findall(r'(MSBNK-[A-Za-z0-9_]+-[A-Za-z0-9]+)', mappingFeatureId)
        accession = accession[0]
        #print(accession, zodiac, row["ZodiacScore"])

        if zodiac > 0.0:
            if not accession in zodiacdict:
                zodiacdict[accession] = []
            temp = {"score":zodiac, "ionform": ionform, "formula": formula}
            zodiacdict[accession].append(temp)

        if rank < rankingthreshold :
            if not accession in siriusdict:
                siriusdict[accession] = []
            temp = {"score":sirius, "ionform": ionform, "formula": formula}
            siriusdict[accession].append(temp)

    for accession, hits in confidencedict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "compound"
            scoreing = "SIRIUS6_Confidence"
            db = "biodatabase"
            score = hit["score"]
            inchi = hit["InChIKey"]
            name = hit["name"]
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, inchi, ionform, name])
            #ddatoolbox.writeLog(logfilename, ["SIRIUS6", id, type, scoreing, db, score, inchi, ionform, name])
            counter += 1
            if counter > rankingthreshold:
                break



    for accession, hits in csidict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "compound"
            scoreing = "SIRIUS6_CSI"
            db = "biodatabase"
            score = hit["score"]
            inchi = hit["InChIKey"]
            name = hit["name"]
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, inchi, ionform, name])
            #ddatoolbox.writeLog(logfilename, ["SIRIUS6", id, type, scoreing, db, score, inchi, ionform, name])
            counter += 1
            if counter > rankingthreshold:
                break




    for accession, hits in zodiacdict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "mf"
            scoreing = "SIRIUS6_Zodiac"
            db = ""
            score = hit["score"]
            formula = hit["formula"]
            name = ""
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, formula, ionform, name])
            #ddatoolbox.writeLog(logfilename, ["SIRIUS6", id, type, scoreing, db, score, formula, ionform, name])

            counter += 1
            if counter > rankingthreshold:
                break

    for accession, hits in siriusdict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "mf"
            scoreing = "SIRIUS6_Sirius"
            db = ""
            score = hit["score"]
            formula = hit["formula"]
            name = ""
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, formula, ionform, name])
            #ddatoolbox.writeLog(logfilename, ["SIRIUS6", id, type, scoreing, db, score, formula, ionform, name])


            counter += 1
            if counter > rankingthreshold:
                break



    with open(outputpath/ "output.txt", mode="a",encoding='UTF-8', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(outputtable)

    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Finished", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "SIRIUS6"])
    ddatoolbox.writeLog(startlogfilename, ["Finished","annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "SIRIUS6"])




def output_sirius5(dataset, path):
    """Parse SIRIUS 5 output directories and write compound/formula annotations to output.txt.

    Iterates over per-accession output directories produced by SIRIUS 5, reads
    ``structure_candidates.tsv`` and ``formula_candidates.tsv`` for each record,
    and extracts Confidence scores, CSI:FingerID scores, Zodiac scores, and Sirius
    scores. All ranked results are appended as tab-separated rows to the dataset's
    ``output.txt`` file.

    Args:
        dataset (str): Dataset identifier string.
        path (dict): Dictionary containing file path information with the following keys:
            ``"dataset"`` (pathlib.Path), ``"software"`` (pathlib.Path),
            ``"logfilename"`` (str), ``"startlogfilename"`` (str),
            ``"output"`` (pathlib.Path), ``"softwareinput"`` (pathlib.Path),
            ``"softwareoutput"`` (pathlib.Path).

    Returns:
        None
    """
    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    softwareinput = path["softwareinput"]
    softwareoutput = path["softwareoutput"]
    #
    #
    # SetUp MSFINDER.INI
    #


    accession_stem_ms2 = dataset+"P"
    accession_stem_ms1 = dataset+"S"


    confidencedict = {}
    csidict = {}

    zodiacdict = {}
    siriusdict = {}

    outputtable = []

    pathlist = list(softwareoutput.glob(r'*'))

    #print(paths, softwareoutput)


    for path in pathlist:
        accession = re.findall(r'(MSBNK-[A-Za-z0-9_]+-[A-Za-z0-9]+)', str(path))
        #
        # Get MassBank record ID
        #
        if len(accession) == 0:
            continue
        accession = accession[0]
        #
        # Compound
        #


        structure_candidates_path = Path(path) / "structure_candidates.tsv"
        #
        # If there is no result
        #
        if not os.path.isfile(structure_candidates_path):
            continue
        #
        # Read structure_candidates.tsv
        #
        rows = ddatoolbox.readtable(structure_candidates_path,\
         ["structureRankPerFormula","molecularFormula", "InChIkey2D","ConfidenceScore", \
         "adduct","CSI:FingerIDScore","name","pubchemids"])

        # If the result is empty
        if len(rows) == 0:
            continue

        for row in rows:
            if row["ConfidenceScore"] == "N/A":
                row["ConfidenceScore"] = float(0.0)
            confidence=float(row["ConfidenceScore"])

            rank = int(row["structureRankPerFormula"])
            csi = float(row["CSI:FingerIDScore"])
            inchikey = row["InChIkey2D"]
            Inchikey1_query = inchikey.split("-")[0]

            ionform=row["adduct"].replace(" ", "")
            formula=row["molecularFormula"]
            name = row["name"]
            pubchem = row["pubchemids"]

            if confidence > 0.0:
                if not accession in confidencedict:
                    confidencedict[accession] = []
                temp = {"score":confidence, "ionform": ionform, "name": name, "formula": formula, "InChIKey": Inchikey1_query}
                confidencedict[accession].append(temp)

            if rank < rankingthreshold :
                if not accession in csidict:
                    csidict[accession] = []
                temp = {"score":csi, "ionform": ionform, "name": name, "formula": formula, "InChIKey": Inchikey1_query}
                csidict[accession].append(temp)

        #
        # mf
        #
        formula_candidates_path = Path(path) / "formula_candidates.tsv"
        #
        # If there is no result
        #
        if not os.path.isfile(formula_candidates_path):
            continue
        #
        # Read formula_candidates.tsv
        #
        rows = ddatoolbox.readtable(formula_candidates_path,\
         ["formulaRank","molecularFormula", "adduct", "SiriusScore","ZodiacScore", "ConfidenceScore"])

        # If the result is empty
        if len(rows) == 0:
            continue

        for row in rows:
            if row["ZodiacScore"] == "N/A":
                row["ZodiacScore"] = float(0.0)
            zodiac=float(row["ZodiacScore"])
            rank = int(row["formulaRank"])
            sirius = float(row["SiriusScore"])
            ionform=row["adduct"].replace(" ", "")
            formula=row["molecularFormula"]

            if zodiac > 0.0:
                if not accession in zodiacdict:
                    zodiacdict[accession] = []
                temp = {"score":zodiac, "ionform": ionform, "formula": formula}
                zodiacdict[accession].append(temp)

            if rank < rankingthreshold :
                if not accession in siriusdict:
                    siriusdict[accession] = []
                temp = {"score":sirius, "ionform": ionform, "formula": formula}
                siriusdict[accession].append(temp)





    for accession, hits in confidencedict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "compound"
            scoreing = "SIRIUS5_Confidence"
            db = "biodatabase"
            score = hit["score"]
            inchi = hit["InChIKey"]
            name = hit["name"]
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, inchi, ionform, name])
            #ddatoolbox.writeLog(logfilename, ["SIRIUS5", id, type, scoreing, db, score, inchi, ionform, name])
            counter += 1
            if counter > rankingthreshold :
                break



    for accession, hits in csidict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "compound"
            scoreing = "SIRIUS5_CSI"
            db = "biodatabase"
            score = hit["score"]
            inchi = hit["InChIKey"]
            name = hit["name"]
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, inchi, ionform, name])
            #ddatoolbox.writeLog(logfilename, ["SIRIUS5", id, type, scoreing, db, score, inchi, ionform, name])
            counter += 1
            if counter > rankingthreshold:
                break




    for accession, hits in zodiacdict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "mf"
            scoreing = "SIRIUS5_Zodiac"
            db = ""
            score = hit["score"]
            formula = hit["formula"]
            name = ""
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, formula, ionform, name])
            #ddatoolbox.writeLog(logfilename, ["SIRIUS5", id, type, scoreing, db, score, formula, ionform, name])

            counter += 1
            if counter > rankingthreshold:
                break

    for accession, hits in siriusdict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "mf"
            scoreing = "SIRIUS5_Sirius"
            db = ""
            score = hit["score"]
            formula = hit["formula"]
            name = ""
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, formula, ionform, name])
            #ddatoolbox.writeLog(logfilename, ["SIRIUS5", id, type, scoreing, db, score, formula, ionform, name])


            counter += 1
            if counter > rankingthreshold:
                break



    with open(outputpath/ "output.txt", mode="a",encoding='UTF-8', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(outputtable)

    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Finished", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "SIRIUS5"])
    ddatoolbox.writeLog(startlogfilename, ["Finished", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "SIRIUS5"])



def setup_sirius(version, dataset, path):
    """Prepare SIRIUS input files (.ms) from MassBank step2 records with correlation filtering.

    Reads MS2 records from the step2 directory, applies correlation-based spectral
    cleaning (``wash_spectra``), recalibrates MS1 m/z values, and writes SIRIUS-format
    ``.ms`` files (one per candidate ion form) to the SIRIUS input directory. Records
    with fewer than 3 peaks after cleaning or without a valid precursor m/z are skipped.
    TOF instruments skip records with precursor m/z above 900.

    Args:
        version (int): SIRIUS version to use (5 or 6); used only for log labeling.
        dataset (str): Dataset identifier string.
        path (dict): Dictionary containing file path information with the following keys:
            ``"dataset"`` (pathlib.Path), ``"software"`` (pathlib.Path),
            ``"logfilename"`` (str), ``"startlogfilename"`` (str),
            ``"output"`` (pathlib.Path), ``"softwareinput"`` (pathlib.Path).

    Returns:
        None
    """
    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    outputpath_SIRIUS = path["softwareinput"]
    #outputpath_ANNOTATION = path["softwareoutput"]
    #
    #
    # SetUp MSFINDER.INI
    #

    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

    #msfinderini = msfinderini.replace("[MS1ttol]", str(settings["ms1_mz_averaged_thres"] ))
    #msfinderini = msfinderini.replace("[MS2ttol]", str(settings["ms2_mz_averaged_thres"] ))

    #accession_stem_ms2 = dataset+"P"
    #accession_stem_ms1 = dataset+"S"

    accession_stem = settings["accession_stem"]

    accession_stem_ms2 = accession_stem+"P"
    accession_stem_ms1 = accession_stem+"S"


    filelist = []
    dt_now = datetime.datetime.now()
    if version == 6:
        siriusversion = "SIRIUS6"
    else:
        siriusversion = "SIRIUS5"

    ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), siriusversion], mode = "a")
    ddatoolbox.writeLog(startlogfilename, ["Started", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), siriusversion], mode = "a")

    files = glob.glob(str(inputpath))

    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)


    for ms2filename in filelist:
        #
        # Loading MassBank recored MS2
        #
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)

        #
        # Getting information
        #
        accession = mb2.get_accession()
        #precmz, intensity = mb2.get_precmz()
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
            continue

        instrument = mb2.get_instrument_type()
        if "TOF" in instrument:
            if float(precmz) > 900.0:
                print(accession, instrument, precmz,  "is skipped")
                continue


        polarity = mb2.get_polarity()
        #retime, sd = mb2.get_rt()
        #retime = retime * 60
        polarity = mb2.get_polarity()
        polsymbol = "+"
        if polarity == "NEGATIVE":
            polarity = "Negative"
            polsymbol = "-"
        if polarity == "POSITIVE":
            polarity = "Positive"
        #
        #PROPOSED_ION_FORM
        #
        Ion_forn_list_temp = mb2.get_proposed_ion_form()
        Ion_forn_list = set()

        if polarity == "Negative":
            Ion_forn_list.add("[M-H]-")
        if polarity == "Positive":
            Ion_forn_list.add("[M+H]+")

        for Ion_forn_temp in Ion_forn_list_temp:
            ionform = Ion_forn_temp.split(" ")[0]
            ionform = "["+ionform+"]"+polsymbol
            ionform = ionform.replace(";", "")
            Ion_forn_list.add(ionform)

        #
        # Spectra cleaning up
        #
        cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)
        if len(cleanedspectra) < 3:
            continue

        #
        #
        # Loading MassBank recored MS1

        mb1 = massbankio.MassBankRecord()
        ms1filename = ms2filename.replace(accession_stem_ms2, accession_stem_ms1)
        ms1filename = ms1filename.replace("step2", "step1")
        #mb1.read_MassBankRecord(ms1filename)
        try:
            mb1.read_MassBankRecord(ms1filename)
        except:
            ms1flag = "failed"
        else:
            ms1flag = "success"

        if ms1flag == "success":

            #
            # Recalibration of ms1 m/z
            #

            settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

            annotationlist = mb1.get_annotation_from_items()
            annoaationlist_corrected = {}
            spec_corrected = {}

            for mzraw, datahash in annotationlist.items():
                mz_corrected = mfsearch.mzcorrection(mzraw, settings["ms1_mz_correction"])
                annoaationlist_corrected[mz_corrected] = datahash
                annoaationlist_corrected[mz_corrected]["mzraw"] = mzraw
                spec_corrected[mz_corrected] = annoaationlist_corrected[mz_corrected]["intensity"]
            mb1.add_annotation_from_items(annoaationlist_corrected, ["m/z(raw)", "ppm(std)", "int", "int(std)", "n", "intlist"])
            mb1.add_spectra(spec_corrected)
            #
            # Spectra cleaning up
            #
            ms1spec = mb1.get_spectra()

        #
        # for SIRIUS5
        #
        #
        counter = 0
        for ionform in Ion_forn_list:
            outputtextSirius = ""
            outputtextSirius += ">compound {}_{}\n".format(accession, counter)
            outputtextSirius += ">ionization {}\n".format(ionform)
            outputtextSirius += ">parentmass {}\n".format(precmz)
            if ms1flag == "success":
                outputtextSirius += ">ms1\n"
                for mz, intens in ms1spec.items():
                    outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            outputtextSirius += ">ms2\n"
            for mz, intens in cleanedspectra.items():
                outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            filenametemp = accession+"_"+str(counter)+".ms"
            with open(outputpath_SIRIUS / filenametemp, mode='w', newline="") as f:
                f.write(outputtextSirius)
            counter += 1




def setup_sirius_nocorrelfilter(version, dataset, path):
    """Prepare SIRIUS input files (.ms) from MassBank step2 records without correlation filtering.

    Similar to ``setup_sirius()``, but uses the raw spectrum (``get_spectra()``) instead
    of applying correlation-based cleaning. Recalibrates MS1 m/z values and writes
    SIRIUS-format ``.ms`` files (one per candidate ion form) to the SIRIUS input directory.
    Records with fewer than 3 peaks or without a valid precursor m/z are skipped.

    Args:
        version (int): SIRIUS version to use (5 or 6); used only for log labeling.
        dataset (str): Dataset identifier string.
        path (dict): Dictionary containing file path information with the following keys:
            ``"dataset"`` (pathlib.Path), ``"software"`` (pathlib.Path),
            ``"logfilename"`` (str), ``"output"`` (pathlib.Path),
            ``"softwareinput"`` (pathlib.Path).

    Returns:
        None
    """
    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    outputpath_SIRIUS = path["softwareinput"]
    #outputpath_ANNOTATION = path["softwareoutput"]
    #
    #
    # SetUp MSFINDER.INI
    #

    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

    #msfinderini = msfinderini.replace("[MS1ttol]", str(settings["ms1_mz_averaged_thres"] ))
    #msfinderini = msfinderini.replace("[MS2ttol]", str(settings["ms2_mz_averaged_thres"] ))

    #accession_stem_ms2 = dataset+"P"
    #accession_stem_ms1 = dataset+"S"

    accession_stem = settings["accession_stem"]

    accession_stem_ms2 = accession_stem+"P"
    accession_stem_ms1 = accession_stem+"S"


    filelist = []
    dt_now = datetime.datetime.now()
    if version == 6:
        siriusversion = "SIRIUS6"
    else:
        siriusversion = "SIRIUS5"

    ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), siriusversion], mode = "a")


    files = glob.glob(str(inputpath))

    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)


    for ms2filename in filelist:
        #
        # Loading MassBank recored MS2
        #
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)

        #
        # Getting information
        #
        accession = mb2.get_accession()
        #precmz, intensity = mb2.get_precmz()
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
            continue

        instrument = mb2.get_instrument_type()
        #if "TOF" in instrument:
        if float(precmz) > 900.0:
            print(accession, instrument, precmz,  "is skipped")
            continue

        #retime, sd = mb2.get_rt()
        #retime = retime * 60
        polarity = mb2.get_polarity()
        polsymbol = "+"
        if polarity == "NEGATIVE":
            polarity = "Negative"
            polsymbol = "-"
        if polarity == "POSITIVE":
            polarity = "Positive"
        #
        #PROPOSED_ION_FORM
        #
        Ion_forn_list_temp = mb2.get_proposed_ion_form()
        Ion_forn_list = set()

        if polarity == "Negative":
            Ion_forn_list.add("[M-H]-")
        if polarity == "Positive":
            Ion_forn_list.add("[M+H]+")

        for Ion_forn_temp in Ion_forn_list_temp:
            ionform = Ion_forn_temp.split(" ")[0]
            ionform = "["+ionform+"]"+polsymbol
            ionform = ionform.replace(";", "")
            Ion_forn_list.add(ionform)

        #
        # Spectra cleaning up
        #
        #cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)
        cleanedspectra = mb2.get_spectra()
        if len(cleanedspectra) < 3:
            continue

        #
        #
        # Loading MassBank recored MS1

        mb1 = massbankio.MassBankRecord()
        ms1filename = ms2filename.replace(accession_stem_ms2, accession_stem_ms1)
        ms1filename = ms1filename.replace("step2", "step1")
        #mb1.read_MassBankRecord(ms1filename)
        try:
            mb1.read_MassBankRecord(ms1filename)
        except:
            ms1flag = "failed"
        else:
            ms1flag = "success"

        if ms1flag == "success":

            #
            # Recalibration of ms1 m/z
            #

            settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

            annotationlist = mb1.get_annotation_from_items()
            annoaationlist_corrected = {}
            spec_corrected = {}

            for mzraw, datahash in annotationlist.items():
                mz_corrected = mfsearch.mzcorrection(mzraw, settings["ms1_mz_correction"])
                annoaationlist_corrected[mz_corrected] = datahash
                annoaationlist_corrected[mz_corrected]["mzraw"] = mzraw
                spec_corrected[mz_corrected] = annoaationlist_corrected[mz_corrected]["intensity"]
            mb1.add_annotation_from_items(annoaationlist_corrected, ["m/z(raw)", "ppm(std)", "int", "int(std)", "n", "intlist"])
            mb1.add_spectra(spec_corrected)
            #
            # Spectra cleaning up
            #
            ms1spec = mb1.get_spectra()

        #
        # for SIRIUS5
        #
        #
        counter = 0
        for ionform in Ion_forn_list:
            outputtextSirius = ""
            outputtextSirius += ">compound {}_{}\n".format(accession, counter)
            outputtextSirius += ">ionization {}\n".format(ionform)
            outputtextSirius += ">parentmass {}\n".format(precmz)
            if ms1flag == "success":
                outputtextSirius += ">ms1\n"
                for mz, intens in ms1spec.items():
                    outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            outputtextSirius += ">ms2\n"
            for mz, intens in cleanedspectra.items():
                outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            filenametemp = accession+"_"+str(counter)+".ms"
            with open(outputpath_SIRIUS / filenametemp, mode='w', newline="") as f:
                f.write(outputtextSirius)
            counter += 1




def setup_sirius_nocal(version, dataset, path):
    """Prepare SIRIUS input files (.ms) using another dataset's proposed ion forms.

    Variant of :func:`setup_sirius` for datasets that lack their own
    ``PROPOSED_ION_FORM`` annotations: correlation-cleaned MS2 spectra are read
    from ``dataset``, but the raw (uncleaned) spectrum is written out and the
    candidate ion-form list is instead read from the corresponding record in a
    fixed reference dataset, ``"HSA001"`` (found by substituting the dataset name
    in the file path). Otherwise follows the same MS1 recalibration and ``.ms``
    file-writing logic as :func:`setup_sirius`.

    Args:
        version (int): SIRIUS version to use (5 or 6); used only for log labeling.
        dataset (str): Dataset identifier string whose MS2 spectra are used.
        path (dict): Dictionary containing file path information (same keys as
            :func:`setup_sirius`).

    Returns:
        None
    """
    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    outputpath_SIRIUS = path["softwareinput"]
    #outputpath_ANNOTATION = path["softwareoutput"]
    #
    #
    # SetUp MSFINDER.INI
    #

    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

    #msfinderini = msfinderini.replace("[MS1ttol]", str(settings["ms1_mz_averaged_thres"] ))
    #msfinderini = msfinderini.replace("[MS2ttol]", str(settings["ms2_mz_averaged_thres"] ))
    accession_stem = settings["accession_stem"]

    accession_stem_ms2 = accession_stem+"P"
    accession_stem_ms1 = accession_stem+"S"
    filelist = []
    dt_now = datetime.datetime.now()
    if version == 6:
        siriusversion = "SIRIUS6"
    else:
        siriusversion = "SIRIUS5"

    ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), siriusversion], mode = "a")


    files = glob.glob(str(inputpath))

    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)


    for ms2filename in filelist:
        #
        # Loading MassBank recored MS2
        #
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)


        #
        # Fetch the corresponding data from HSA001
        #
        #print(ms2filename)

        ms2filename = ms2filename.replace(dataset, "HSA001")
        #print(dataset, ms2filename)

        mb21 = massbankio.MassBankRecord()
        mb21.read_MassBankRecord(ms2filename)


        #
        # Getting information
        #
        accession = mb2.get_accession()
        #precmz, intensity = mb2.get_precmz()
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
            continue

        #retime, sd = mb2.get_rt()
        #retime = retime * 60
        polarity = mb2.get_polarity()
        polsymbol = "+"
        if polarity == "NEGATIVE":
            polarity = "Negative"
            polsymbol = "-"
        if polarity == "POSITIVE":
            polarity = "Positive"
        #
        #PROPOSED_ION_FORM
        #
        # Adduct data comes from HSA001
        Ion_forn_list_temp = mb21.get_proposed_ion_form()
        Ion_forn_list = set()

        if polarity == "Negative":
            Ion_forn_list.add("[M-H]-")
        if polarity == "Positive":
            Ion_forn_list.add("[M+H]+")

        for Ion_forn_temp in Ion_forn_list_temp:
            ionform = Ion_forn_temp.split(" ")[0]
            ionform = "["+ionform+"]"+polsymbol
            ionform = ionform.replace(";", "")
            Ion_forn_list.add(ionform)

        #
        # Spectra cleaning up
        #
        cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)
        if len(cleanedspectra) < 3:
            continue

        cleanedspectra = mb2.get_spectra()



        #
        #
        # Loading MassBank recored MS1

        mb1 = massbankio.MassBankRecord()
        ms1filename = ms2filename.replace(accession_stem_ms2, accession_stem_ms1)
        ms1filename = ms1filename.replace("step2", "step1")
        #mb1.read_MassBankRecord(ms1filename)
        try:
            mb1.read_MassBankRecord(ms1filename)
        except:
            ms1flag = "failed"
        else:
            ms1flag = "success"

        if ms1flag == "success":

            #
            # Recalibration of ms1 m/z
            #

            settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

            annotationlist = mb1.get_annotation_from_items()
            annoaationlist_corrected = {}
            spec_corrected = {}

            for mzraw, datahash in annotationlist.items():
                mz_corrected = mfsearch.mzcorrection(mzraw, settings["ms1_mz_correction"])
                annoaationlist_corrected[mz_corrected] = datahash
                annoaationlist_corrected[mz_corrected]["mzraw"] = mzraw
                spec_corrected[mz_corrected] = annoaationlist_corrected[mz_corrected]["intensity"]
            mb1.add_annotation_from_items(annoaationlist_corrected, ["m/z(raw)", "ppm(std)", "int", "int(std)", "n", "intlist"])
            mb1.add_spectra(spec_corrected)
            #
            # Spectra cleaning up
            #
            ms1spec = mb1.get_spectra()

        #
        # for SIRIUS5
        #
        #
        counter = 0
        for ionform in Ion_forn_list:
            outputtextSirius = ""
            outputtextSirius += ">compound {}_{}\n".format(accession, counter)
            outputtextSirius += ">ionization {}\n".format(ionform)
            outputtextSirius += ">parentmass {}\n".format(precmz)
            if ms1flag == "success":
                outputtextSirius += ">ms1\n"
                for mz, intens in ms1spec.items():
                    outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            outputtextSirius += ">ms2\n"
            for mz, intens in cleanedspectra.items():
                outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            filenametemp = accession+"_"+str(counter)+".ms"
            with open(outputpath_SIRIUS / filenametemp, mode='w', newline="") as f:
                f.write(outputtextSirius)
            counter += 1




def setup_sirius_raw(version, dataset, path):
    """Prepare SIRIUS input files (.ms) using one representative raw (un-averaged) spectrum.

    Experimental variant of :func:`setup_sirius`: instead of using the averaged
    MS1/MS2 spectra stored in the MassBank record, this picks one representative
    raw scan (the base-peak-intensity-ranked scan at the 20th percentile,
    ``rawspec_middle_key``) referenced in the record's ``COMMENT: RAW_SPEC``
    field, loads the corresponding raw Hive (``.hmd``) file directly via the
    ``hive_data_access`` package, and writes its raw MS1/MS2 peaks to the
    ``.ms`` file instead of the averaged spectra.

    Args:
        version (int): SIRIUS version to use (5 or 6); used only for log labeling.
        dataset (str): Dataset identifier string.
        path (dict): Dictionary containing file path information (same keys as
            :func:`setup_sirius`).

    Returns:
        None

    Note:
        Requires the external ``hive_data_access`` package (not part of
        ``ddatoolbox``) to be importable. Contains dataset-name substitution
        rules (``"HSA2"->"HSA0"``, etc.) hard-coded for a specific raw-data
        folder layout.
    """
    #from hivetools import HiveFile, ValueRange
    from hive_data_access import HiveFile, ValueRange
    from . import ddatoolbox

    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    outputpath_SIRIUS = path["softwareinput"]

    #outputpath_ANNOTATION = path["softwareoutput"]
    #
    #
    # SetUp MSFINDER.INI
    #

    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

    #msfinderini = msfinderini.replace("[MS1ttol]", str(settings["ms1_mz_averaged_thres"] ))
    #msfinderini = msfinderini.replace("[MS2ttol]", str(settings["ms2_mz_averaged_thres"] ))
    accession_stem = settings["accession_stem"]

    accession_stem_ms2 = accession_stem+"P"
    accession_stem_ms1 = accession_stem+"S"


    filelist = []
    dt_now = datetime.datetime.now()
    if version == 6:
        siriusversion = "SIRIUS6"
    else:
        siriusversion = "SIRIUS5"

    #ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), siriusversion], mode = "a")


    files = glob.glob(str(inputpath))

    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)


    for ms2filename in filelist:

        #
        # Loading MassBank recored MS2
        #
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)

        #
        # Getting information
        #
        accession = mb2.get_accession()
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
            continue

        instrument = mb2.get_instrument_type()
        #if "TOF" in instrument:
        if float(precmz) > 900.0:
            print(accession, instrument, precmz,  "is skipped")
            continue
        #precmz, intensity = mb2.get_precmz()
        #retime, sd = mb2.get_rt()
        #retime = retime * 60
        polarity = mb2.get_polarity()
        polsymbol = "+"
        if polarity == "NEGATIVE":
            polarity = "Negative"
            polsymbol = "-"
        if polarity == "POSITIVE":
            polarity = "Positive"
        #
        #PROPOSED_ION_FORM
        #
        Ion_forn_list_temp = mb2.get_proposed_ion_form()
        Ion_forn_list = set()

        if polarity == "Negative":
            Ion_forn_list.add("[M-H]-")
        if polarity == "Positive":
            Ion_forn_list.add("[M+H]+")

        for Ion_forn_temp in Ion_forn_list_temp:
            ionform = Ion_forn_temp.split(" ")[0]
            ionform = "["+ionform+"]"+polsymbol
            ionform = ionform.replace(";", "")
            Ion_forn_list.add(ionform)

        #
        # Spectra cleaning up
        #
        cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)
        if len(cleanedspectra) < 3:
            continue



        rawdata = eval(mb2.get_anything("COMMENT:", "RAW_SPEC"))
        rawspeckeys = list(rawdata.keys())
        rawspeckeys_sorted = sorted(rawspeckeys, key = lambda x: rawdata[x]['bp'])
        rawspec_middle_key = rawspeckeys_sorted[int(len(rawspeckeys_sorted)/5)]

        print(accession, rawspec_middle_key)
        filename_scan = rawspec_middle_key.split("_")
        scan = int(filename_scan[-1])
        filename = "_".join(filename_scan[0:-1])+".hmd"
        #print(filename_scan, filename, scan)
        rawdatapath = datasetpath / "hmd" / filename

        rawdatapath = str(rawdatapath)
        rawdatapath = rawdatapath.replace("HSA2","HSA0")
        rawdatapath = rawdatapath.replace("HSA3","HSA1")
        rawdatapath = rawdatapath.replace("HSA4","HSA0")
        rawdatapath = rawdatapath.replace("HSA6","HSA1")
        #print(rawdatapath)

        with HiveFile(str(rawdatapath)) as hive_file:

            measurement = hive_file.get_measurement(0)
            qual_analysis = measurement.get_qual_analysis(0)
            spectrum = qual_analysis.get_spectrum(scan)

            precursormz = spectrum.acquisition_condition.precursor_value
            ms_level = str(spectrum.acquisition_condition.ms_level)
            retentionTime = str(spectrum.retention_time)
            precursorIndex = str(spectrum.precursor_index)


            msmsdict, basePeakIntensity, norm = ddatoolbox.generateSpecHash(spectrum, threshold = 0)



            spectrum = qual_analysis.get_spectrum(int(precursorIndex))

            ms1retentionTime = str(spectrum.retention_time)
            ms1dict, ms1basePeakIntensity, ms1norm = ddatoolbox.generateSpecHash(spectrum, threshold = 0)

        #print(msmsdict)
        #print(ms1dict)

        #
        # Spectra cleaning up
        #
        #ms1spec = mb1.get_spectra()
        ms1spec = ms1dict
        precmz = precursormz

        #
        # for SIRIUS5
        #
        #

        cleanedspectra = msmsdict
        #cleanedspectra = mb2.get_spectra()

        if len(cleanedspectra) < 3:
            continue

        counter = 0
        Ion_forn_list_sorted = sorted(Ion_forn_list)
        for ionform in Ion_forn_list_sorted:
            outputtextSirius = ""
            outputtextSirius += ">compound {}_{}\n".format(accession, counter)
            outputtextSirius += ">ionization {}\n".format(ionform)
            outputtextSirius += ">parentmass {}\n".format(precmz)
            outputtextSirius += ">ms1\n"
            for mz, intens in ms1spec.items():
                outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            outputtextSirius += ">ms2\n"
            for mz, intens in cleanedspectra.items():
                outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            filenametemp = accession+"_"+str(counter)+".ms"
            with open(outputpath_SIRIUS / filenametemp, mode='w', newline="") as f:
                f.write(outputtextSirius)
            counter += 1




def setup_sirius_raw2(version, dataset_original, path):
    """Prepare SIRIUS input files (.ms) from a raw scan, resolved via an explicit source dataset.

    Variant of :func:`setup_sirius_raw` that resolves the raw Hive (``.hmd``)
    file path directly as ``../data/<dataset_original>/hmd/<filename>`` (no
    dataset-name substitution rules) and skips records whose raw file cannot be
    loaded (wrapped in a try/except) rather than failing outright. No TOF
    precursor-m/z upper-limit filtering is applied.

    Args:
        version (int): SIRIUS version to use (5 or 6); used only for log labeling.
        dataset_original (str): Dataset identifier whose raw ``hmd/`` folder
            contains the referenced raw scan file.
        path (dict): Dictionary containing file path information (same keys as
            :func:`setup_sirius`).

    Returns:
        None

    Note:
        Requires the external ``hive_data_access`` package to be importable.
    """
    #from hivetools import HiveFile, ValueRange
    from hive_data_access import HiveFile, ValueRange
    from . import ddatoolbox

    #datasetname = "HSA{:0>3}"
    #dataset = datasetname.format(i*10 + 101)
    #dataset_original = datasetname.format(i)


    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    outputpath_SIRIUS = path["softwareinput"]

    #outputpath_ANNOTATION = path["softwareoutput"]
    #
    #
    # SetUp MSFINDER.INI
    #

    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

    #msfinderini = msfinderini.replace("[MS1ttol]", str(settings["ms1_mz_averaged_thres"] ))
    #msfinderini = msfinderini.replace("[MS2ttol]", str(settings["ms2_mz_averaged_thres"] ))
    accession_stem = settings["accession_stem"]

    accession_stem_ms2 = accession_stem+"P"
    accession_stem_ms1 = accession_stem+"S"


    filelist = []
    dt_now = datetime.datetime.now()
    if version == 6:
        siriusversion = "SIRIUS6"
    else:
        siriusversion = "SIRIUS5"

    #ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), siriusversion], mode = "a")


    files = glob.glob(str(inputpath))

    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)


    for ms2filename in filelist:

        #
        # Loading MassBank recored MS2
        #
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)

        #
        # Getting information
        #
        accession = mb2.get_accession()
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
            continue
        #precmz, intensity = mb2.get_precmz()
        #retime, sd = mb2.get_rt()
        #retime = retime * 60
        polarity = mb2.get_polarity()
        polsymbol = "+"
        if polarity == "NEGATIVE":
            polarity = "Negative"
            polsymbol = "-"
        if polarity == "POSITIVE":
            polarity = "Positive"
        #
        #PROPOSED_ION_FORM
        #
        Ion_forn_list_temp = mb2.get_proposed_ion_form()
        Ion_forn_list = set()

        if polarity == "Negative":
            Ion_forn_list.add("[M-H]-")
        if polarity == "Positive":
            Ion_forn_list.add("[M+H]+")

        for Ion_forn_temp in Ion_forn_list_temp:
            ionform = Ion_forn_temp.split(" ")[0]
            ionform = "["+ionform+"]"+polsymbol
            ionform = ionform.replace(";", "")
            Ion_forn_list.add(ionform)

        #
        # Spectra cleaning up
        #
        cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)
        if len(cleanedspectra) < 3:
            continue



        rawdata = eval(mb2.get_anything("COMMENT:", "RAW_SPEC"))
        rawspeckeys = list(rawdata.keys())
        rawspeckeys_sorted = sorted(rawspeckeys, key = lambda x: rawdata[x]['bp'])
        rawspec_middle_key = rawspeckeys_sorted[int(len(rawspeckeys_sorted)/5)]

        print(accession, rawspec_middle_key)
        filename_scan = rawspec_middle_key.split("_")
        scan = int(filename_scan[-1])
        filename = "_".join(filename_scan[0:-1])+".hmd"
        #print(filename_scan, filename, scan)
        #rawdatapath = datasetpath / "hmd" / filename
        rawdatapath = "../data/"+dataset_original+"/hmd/"+filename
        print(rawdatapath)


        try:

            with HiveFile(str(rawdatapath)) as hive_file:

                measurement = hive_file.get_measurement(0)
                qual_analysis = measurement.get_qual_analysis(0)
                spectrum = qual_analysis.get_spectrum(scan)

                precursormz = spectrum.acquisition_condition.precursor_value
                ms_level = str(spectrum.acquisition_condition.ms_level)
                retentionTime = str(spectrum.retention_time)
                precursorIndex = str(spectrum.precursor_index)

                msmsdict, basePeakIntensity, norm = ddatoolbox.generateSpecHash(spectrum, threshold = 0)

                spectrum = qual_analysis.get_spectrum(int(precursorIndex))

                ms1retentionTime = str(spectrum.retention_time)
                ms1dict, ms1basePeakIntensity, ms1norm = ddatoolbox.generateSpecHash(spectrum, threshold = 0)
        except:
            continue

        #print(msmsdict)
        #print(ms1dict)

        #
        # Spectra cleaning up
        #
        #ms1spec = mb1.get_spectra()
        ms1spec = ms1dict
        precmz = precursormz

        #
        # for SIRIUS5
        #
        #

        cleanedspectra = msmsdict
        #cleanedspectra = mb2.get_spectra()

        if len(cleanedspectra) < 3:
            continue

        counter = 0
        Ion_forn_list_sorted = sorted(Ion_forn_list)
        for ionform in Ion_forn_list_sorted:
            outputtextSirius = ""
            outputtextSirius += ">compound {}_{}\n".format(accession, counter)
            outputtextSirius += ">ionization {}\n".format(ionform)
            outputtextSirius += ">parentmass {}\n".format(precmz)
            outputtextSirius += ">ms1\n"
            for mz, intens in ms1spec.items():
                outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            outputtextSirius += ">ms2\n"
            for mz, intens in cleanedspectra.items():
                outputtextSirius += "{:.6f} {}\n".format(mz, intens )
            filenametemp = accession+"_"+str(counter)+".ms"
            with open(outputpath_SIRIUS / filenametemp, mode='w', newline="") as f:
                f.write(outputtextSirius)
            counter += 1


def search_msfinder(dataset, path):
    """Run MS-FINDER over prepared input files and parse its output into ``output.txt``.

    Builds ``MSFINDER.INI`` from the bundled template (substituting the dataset's
    MS1/MS2 averaged m/z thresholds), invokes the MS-FINDER CLI's ``predict`` and
    ``mssearch`` commands over the tool's input/output folders, then parses the
    resulting ``Struc*`` (structure), ``Formula*``, and ``MsSearch*`` result files
    into compound/class/formula/spectral-similarity annotation rows, appended to
    the dataset's ``output.txt``.

    Args:
        dataset (str): Dataset identifier, used to build the accession-stem prefix
            for logging.
        path (dict): Dictionary containing file path information with keys
            ``"dataset"``, ``"software"``, ``"logfilename"``,
            ``"startlogfilename"``, ``"output"``, ``"softwareinput"``
            (MS-FINDER's ``.mat`` input folder), and ``"softwareoutput"``
            (MS-FINDER's result folder).

    Returns:
        None

    Note:
        The block that generates the ``.mat`` peak-list input files (from the
        dataset's step2 MassBank records) is currently disabled (kept as
        commented-out dead code); this function only configures
        ``MSFINDER.INI`` and invokes MS-FINDER, so the ``.mat`` files under
        ``path["softwareinput"]`` must already exist (e.g. from a previous run
        or an external step) before calling this function.
    """
    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    outputpath_MSFINDER = path["softwareinput"]
    outputpath_ANNOTATION = path["softwareoutput"]
    #
    #
    # SetUp MSFINDER.INI
    #

    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")
    with open(softwareroot/ 'MSFINDER_template.txt', 'r+') as f:
        msfinderini = f.read()
    msfinderini = msfinderini.replace("[MS1ttol]", str(settings["ms1_mz_averaged_thres"] ))
    msfinderini = msfinderini.replace("[MS2ttol]", str(settings["ms2_mz_averaged_thres"] ))
    #print(msfinderini)

    pathtomsfinderini = outputpath_MSFINDER / "MSFINDER.INI"

    with open(pathtomsfinderini, 'w') as f:
        f.write(msfinderini)
    #
    accession_stem_ms2 = dataset+"P"
    accession_stem_ms1 = dataset+"S"
    filelist = []
    counterall = 0

    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "MS-FINDER"], mode = "a")
    ddatoolbox.writeLog(startlogfilename, ["Started", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "MS-FINDER"], mode = "a")

    files = glob.glob(str(inputpath))


    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)


    for ms2filename in filelist:
        #
        # Loading MassBank recored MS2
        #
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)
        #
        # Getting information
        #
        accession = mb2.get_accession()
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
            continue

        #retime, sd = mb2.get_rt()
        #retime = retime * 60
        retime = 1.23
        polarity = mb2.get_polarity()
        #print(polarity)
        polsymbol = "+"
        if polarity == "NEGATIVE":
            polarity = "Negative"
            polsymbol = "-"
        if polarity == "POSITIVE":
            polarity = "Positive"
        #
        #PROPOSED_ION_FORM
        #
        Ion_forn_list_temp = mb2.get_proposed_ion_form()
        Ion_forn_list = set()

        if polarity == "Negative":
            Ion_forn_list.add("[M-H]-")
        if polarity == "Positive":
            Ion_forn_list.add("[M+H]+")

        for Ion_forn_temp in Ion_forn_list_temp:
            ionform = Ion_forn_temp.split("; ")[0]
            ionform = "["+ionform+"]"+polsymbol
            Ion_forn_list.add(ionform)

        #
        # Spectra cleaning up by correlation coefficient to precursor
        #
        cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)

        #annotationlist = mb2.get_annotation_from_items()
        #for mz, temp in annotationlist.items():
        #    intensity = temp['intensity']
        #    correlation = temp['correlation']
        #    number = temp['number']
        #    annotationlist[mz] = {"intensity": intensity, "correlation": correlation, "number": number}
        #cleanedspectra = {}
        #for mzraw, datahash in annotationlist.items():
        #    if datahash["correlation"] > production_correlation_threshold:
        #        cleanedspectra[mzraw] = datahash["intensity"]
        if len(cleanedspectra) < 3:
            continue

        #
        #
        # Loading MassBank recored MS1

        mb1 = massbankio.MassBankRecord()
        ms1filename = ms2filename.replace(accession_stem_ms2, accession_stem_ms1)
        ms1filename = ms1filename.replace("step2", "step1")
        try:
            mb1.read_MassBankRecord(ms1filename)
        except:
            ms1flag = "failed"
        else:
            ms1flag = "success"

        if ms1flag == "success":

            #
            # Recalibration of ms1 m/z
            #
            settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")

            annotationlist = mb1.get_annotation_from_items()
            annoaationlist_corrected = {}
            spec_corrected = {}

            for mzraw, datahash in annotationlist.items():
                mz_corrected = mfsearch.mzcorrection(mzraw, settings["ms1_mz_correction"])
                annoaationlist_corrected[mz_corrected] = datahash
                annoaationlist_corrected[mz_corrected]["mzraw"] = mzraw
                spec_corrected[mz_corrected] = annoaationlist_corrected[mz_corrected]["intensity"]
            mb1.add_annotation_from_items(annoaationlist_corrected, ["m/z(raw)", "ppm(std)", "int", "int(std)", "n", "intlist"])
            mb1.add_spectra(spec_corrected)
            #
            # Spectra cleaning up
            #
            ms1spec = mb1.get_spectra()


        counterall+=1

        #
        # for MS-FINDER
        #
        #
        counter = 0
        print(accession, Ion_forn_list)
        for ionform in Ion_forn_list:
            outputtextMSFINDER = ""
            outputtextMSFINDER += "NAME: {}\n".format(accession)
            outputtextMSFINDER += "RETENTIONTIME: {}\n".format(retime)
            outputtextMSFINDER += "PRECURSORMZ: {}\n".format(precmz)
            outputtextMSFINDER += "PRECURSORTYPE: {}\n".format(ionform)
            outputtextMSFINDER += "IONMODE: {}\n".format(polarity)
            if ms1flag == "success":
                outputtextMSFINDER += "MSTYPE: MS1\n"
                outputtextMSFINDER += "Num Peaks: {}\n".format(len(ms1spec))
                for mz, intens in ms1spec.items():
                    outputtextMSFINDER += "{:.6f} {}\n".format(mz, intens)
            outputtextMSFINDER += "MSTYPE: MS2\n"
            outputtextMSFINDER += "Num Peaks: {}\n".format(len(cleanedspectra))
            for mz, intens in cleanedspectra.items():
                outputtextMSFINDER += "{:.6f} {}\n".format(mz, intens)
            tempfilename = accession+"_"+str(counter)+".mat"
            with open(outputpath_MSFINDER / tempfilename, mode='w', newline="") as f:
                f.write(outputtextMSFINDER)
            counter += 1

    
    execmsfinder = paths.getsettings("SOFTWAREPATH", "execmsfinder")
    execmsfinder = Path(execmsfinder).resolve()
    command = r"{exec:s} predict -i {input:s} -o {output:s} -m {ini:s}"
    command = command.format( exec = str(execmsfinder),
    input = str(outputpath_MSFINDER),
    output= str(outputpath_ANNOTATION),
    ini = str(pathtomsfinderini)
    )
    print(command)
    subprocess.run(command, shell=True, cwd=str(outputpath_ANNOTATION))


    #command = execmsfinder + " mssearch -i {input:s} -o {output:s} -m {ini:s}"
    command = r"{exec:s} mssearch -i {input:s} -o {output:s} -m {ini:s}"
    command = command.format(exec = str(execmsfinder),
    input = str(outputpath_MSFINDER),
    output= str(outputpath_ANNOTATION),
    ini = str(pathtomsfinderini)
    )
    print(command)
    subprocess.run(command, shell=True, cwd=str(outputpath_ANNOTATION))



    #file = glob.glob(Path_list.ms_finder_result_upload_path()+"\\*")
    file = glob.glob(str(outputpath_ANNOTATION / "*"))
    result_integral_text_box = []
    for pathdata in file:
        if "Struc" in pathdata:


            with open(pathdata, "r", encoding="utf-8") as f:
                result_data = f.readlines()

                thisline = []

                for text_data in result_data[1:]:
                    converted_list_data = text_data.split("\t")
                    id = converted_list_data[2]
                    type = "compound"
                    scoreing = "MS-FINDER"
                    db = "MS-FINDER_DB"
                    score = converted_list_data[8]
                    inchi = converted_list_data[12]
                    name = converted_list_data[7]
                    ionform = converted_list_data[6]
                    ontology = converted_list_data[11]

                    result_integral_text_box.append([id, type, scoreing, db, score, inchi, ionform, name])
                    #ddatoolbox.writeLog(logfilename, ["MS-FINDER", id, type, scoreing, db, score, inchi, ionform, name])
                    if not ontology == "":
                        type = "class"
                        result_integral_text_box.append([id, type, scoreing, db, score, inchi, ionform, ontology])
                        #ddatoolbox.writeLog(logfilename, ["MS-FINDER", id, type, scoreing, db, score, inchi, ionform, ontology])


        if "Formula" in pathdata:

            with open(pathdata, "r", encoding="utf-8") as f:
                result_data = f.readlines()

                thisline = []

                for text_data in result_data[1:]:
                    converted_list_data = text_data.split("\t")
                    id = converted_list_data[2]
                    type = "mf"
                    scoreing = "MS-FINDER"
                    db = ""
                    score = converted_list_data[10]
                    mf = converted_list_data[7]
                    diff = converted_list_data[9]
                    ionform = converted_list_data[6]
                    result_integral_text_box.append([id, type, scoreing, db, score, mf, ionform, diff])

        if "MsSearch" in pathdata:


            with open(pathdata, "r", encoding="utf-8") as f:
                result_data = f.readlines()
                result_data

                thisline = []

                for text_data in result_data[1:]:
                    converted_list_data = text_data.split("\t")
                    id = converted_list_data[1]
                    #id = id.replace(r"_[0-9]+$", "")
                    id = re.sub(r"_[0-9]+$", "", id)
                    type = "similarity"
                    scoreing = "MS-FINDER"
                    db = "MS-FINDER_DB"
                    score = converted_list_data[9].rstrip()
                    inchi = converted_list_data[7]
                    #name = converted_list_data[7]
                    ionform = converted_list_data[5]

                    result_integral_text_box.append([id, type, scoreing, db, score, inchi, ionform, ""])




    with open(outputpath/ "output.txt", mode="a",encoding='UTF-8', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(result_integral_text_box)

    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Finished", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "MS-FINDER"])
    ddatoolbox.writeLog(startlogfilename, ["Finished", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "MS-FINDER"])


def search_metfrag(dataset, path, database = "PubChemLite"):
    """Run MetFrag against a local structure database for each formula-annotated record.

    Loads the requested local candidate database (either the MS-FINDER structure
    DB, augmented with a difference-library user DB, or PubChemLite), then for
    every step2 MassBank record with at least one supported precursor ion form
    and formula (from ``COMMENT: PROPOSED_MF_FROM_PRECMZ``), builds a per-formula
    candidate subset (filtered by exact formula match, or by m/z tolerance for
    the two special-cased single-adduct ion forms), writes it to
    ``candidates.csv``, fills in the MetFrag parameter template, and invokes the
    MetFrag command-line jar via subprocess. After all records are processed,
    parses every MetFrag result CSV in the output folder and appends ranked
    compound annotation rows to the dataset's ``output.txt``.

    Args:
        dataset (str): Dataset identifier, used to build the accession-stem
            prefix for locating step2 record files.
        path (dict): Dictionary containing file path information with keys
            ``"dataset"``, ``"software"``, ``"logfilename"``,
            ``"startlogfilename"``, ``"output"``, ``"softwareinput"``, and
            ``"softwareoutput"``.
        database (str, optional): Which local structure database to search:
            ``"msfinder"`` (MS-FINDER structure DB), ``"PubChemLite"``, or any
            other value (generic ``"MetFrag"`` labeling). Defaults to
            ``"PubChemLite"``.

    Returns:
        None
    """
    #
    # database = "msfinder"

    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    softwareinput = path["softwareinput"]
    softwareoutput = path["softwareoutput"]
    #
    accession_stem_ms2 = dataset+"P"
    accession_stem_ms1 = dataset+"S"

    if database == "msfinder":
        softwarename = "MetFragMSFINDER"
    elif database == "PubChemLite":
        softwarename = "MetFragPubChemLite"
    else:
        softwarename = "MetFrag"
    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), softwarename], mode = "a")
    ddatoolbox.writeLog(startlogfilename, ["Started", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), softwarename], mode = "a")
    #
    #
    #


    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")
    """
    MetFrag_pubchemlite_db =  ../software/PubChemLite_exposomics_20241227.csv
    MetFrag_pubchem_template =  ../software/MetFragParameterPubChemLite_template.txt
    MetFrag_msfinder_db =  ../software/MsfinderStructureDB-VS15_forMetFrag.csv
    MetFrag_msfinder_template =  ../software/MetFragParameterMSFINDER_template.txt
    MetFrag_msfinder_userdb =  ../software/difflibrary_top2_MetFrag.txt
    """
    #LocalDatabasePath = softwareroot/"kegg_2017-07-23.csv"
    if database == "msfinder":
        #MetFragParameter_template = "MetFragParameterMSFINDER_template.txt"
        MetFragParameter_template = paths.getsettings("SOFTWAREPATH", "MetFrag_msfinder_template")
        LocalDatabasePath = paths.getsettings("SOFTWAREPATH", "MetFrag_msfinder_db") #] softwareroot/"MsfinderStructureDB-VS15_forMetFrag.csv"
        useddb = "MS-FINDER_DB"
        headerlist = ["Identifier", "InChI", "MonoisotopicMass", "MolecularFormula", "InChIKey1", "InChIKey2", "SMILES", "Name", "InChIKey3"]

        print(LocalDatabasePath)
        csv_db = ddatoolbox.readsmalltable(LocalDatabasePath, keywords = headerlist, filetype = "csv")

        #
        # added 20251207
        #
        filenametemp = paths.getsettings("SOFTWAREPATH", "MetFrag_msfinder_userdb")
        filenametemp = str(Path(filenametemp).resolve())
        print(filenametemp)
        csv_db_diff = ddatoolbox.readsmalltable(filenametemp, keywords = headerlist, filetype = "tab")
        csv_db.extend(csv_db_diff)


    elif database == "PubChemLite":
        #
        # Load compound db
        #
        LocalDatabasePath = paths.getsettings("SOFTWAREPATH", "MetFrag_pubchemlite_db") #softwareroot/"PubChemLite_exposomics_20241227.csv"
        useddb = "PubChemLite"
        MetFragParameter_template = paths.getsettings("SOFTWAREPATH", "MetFrag_pubchemlite_template") #"MetFragParameterPubChemLite_template.txt"
        headerlist = ["Identifier","FirstBlock","PubMed_Count","Patent_Count","Related_CIDs","Synonym","MolecularFormula","SMILES","InChI","InChIKey","MonoisotopicMass","CompoundName","AnnoTypeCount","Identification"]

        print(LocalDatabasePath)
        csv_db = ddatoolbox.readsmalltable(LocalDatabasePath, keywords = headerlist, filetype = "csv")
    #csv_db = sorted(csv_db, key = lambda x: x["Exact mass"])

    #"Identifier", "InChI", "MonoisotopicMass", "MolecularFormula", "InChIKey1", "InChIKey2", "SMILES", "Name", "InChIKey3"

    LocalDatabasePath = str(Path(LocalDatabasePath).resolve())
    MetFragParameter_template= str(Path(MetFragParameter_template).resolve())

    filelist = []
    counterall = 0

    ionformdict ={}

    ionformnumbers = {
    "[M+H]+":"1",
    "[M+NH4]+":"18",
    "[M+Na]+":"23",
    "[M+K]+":"39",
    "[M+CH3OH+H]+":"33",
    "[M+ACN+H]+":"42",
    "[M+ACN+Na]+":"64",
    "[M+2ACN+H]+":"83",
    "[M-H]-":"-1",
    "[M+Cl]-":"35",
    "[M+FA-H]-":"45",
    "[M+Hac-H]-":"59"
    }


    files = glob.glob(str(inputpath))

    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)

    for ms2filename in filelist:
        #if not "HSA057P0429007" in ms2filename:
        #    continue
        #
        # Loading MassBank recored MS2
        #
        mb2 = ddatoolbox.massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)
        #
        # Getting information
        #
        accession = mb2.get_accession()
        #precmz, intensity = mb2.get_precmz()
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
            continue

        ms1thres = float(mb2.get_anything("COMMENT:","MS1_M/Z_AVERAGED_THRES"))
        ms2thres = float(mb2.get_anything("COMMENT:","MS2_M/Z_AVERAGED_THRES"))

        polarity = mb2.get_polarity()

        if polarity == "NEGATIVE":
            polsymbol = "-"
            positivecheck = "False"
        else:
            polsymbol = "+"
            positivecheck = "True"
        #
        #PROPOSED_ION_FORM
        #
        Ion_forn_list_temp = mb2.get_anything("COMMENT:", "PROPOSED_MF_FROM_PRECMZ")
        Ion_forn_list = {}


        for Ion_forn_temp in Ion_forn_list_temp:
            mf = Ion_forn_temp.split(" ")[0]
            ionform = Ion_forn_temp.split(" ")[2]
            ionform = "["+ionform+"]"+polsymbol
            ionform = ionform.replace(";", "")


            if ionform in ionformnumbers:
                Ion_forn_list[mf] = ionform
        if len(Ion_forn_list) == 0:
            print("No mf", accession)
            continue
        #
        # Spectra cleaning up
        #
        cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)


        if len(cleanedspectra) < 3:
            continue

        #
        #
        # MetFrag tempspec creation
        #
        outputtext = ""
        for mz, intens in cleanedspectra.items():
            outputtext += "{}\t{}\n".format(mz, intens)

        with open(softwareinput / "tempspec.txt", 'w', encoding="utf-8", newline="") as f:
            f.write(outputtext)
        #
        # for MetFrag
        #
        #
        counter = 0

        print(accession, Ion_forn_list)
        for formula, ionform in Ion_forn_list.items():

            candidates_lst = [headerlist]

            for datatemp in csv_db:
                MolecularFormula = datatemp["MolecularFormula"]
                if MolecularFormula == formula:

                    templist = [datatemp[x] for x in headerlist]
                    candidates_lst.append(templist)
                elif ionform == "[M-H]-":
                    compoundmw = precmz + 1.0072764
                    absolutethres = compoundmw * 0.000001 * ms1thres * 2
                    compoundmwub = compoundmw + absolutethres
                    compoundmwlb = compoundmw - absolutethres
                    em = float(datatemp["MonoisotopicMass"])
                    if compoundmwlb < em < compoundmwub:
                        templist = [datatemp[x] for x in headerlist]
                        candidates_lst.append(templist)
                elif ionform == "[M+H]+":
                    compoundmw = precmz - 1.0072764
                    absolutethres = compoundmw * 0.000001 * ms1thres * 2
                    compoundmwub = compoundmw + absolutethres
                    compoundmwlb = compoundmw - absolutethres
                    em = float(datatemp["MonoisotopicMass"])
                    if compoundmwlb < em < compoundmwub:
                        templist = [datatemp[x] for x in headerlist]
                        candidates_lst.append(templist)


            filenamecandidates = softwareinput/'candidates.csv'

            with open(filenamecandidates, mode="w", encoding='UTF-8', newline="") as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerows(candidates_lst)

            ionformdict[accession +"_" +str(counter)] = ionform
            #
            # for MetFrag
            #
            #
            #with open(softwareroot / MetFragParameter_template, 'r+', encoding="utf-8") as f:
            with open(MetFragParameter_template, 'r+', encoding="utf-8") as f:
                parameters = f.read()
            parameters = parameters.replace("[MS2ttol]", str(ms2thres))
            parameters = parameters.replace("[inputpath]", str(softwareinput))
            parameters = parameters.replace("[outputpath]", str(softwareoutput))
            parameters = parameters.replace("[outputname]", accession +"_" +str(counter))
            parameters = parameters.replace("[formula]",formula)


            parameters = parameters.replace("[precmz]", str(precmz))
            parameters = parameters.replace("[ionform]",str(ionformnumbers[ionform]))
            parameters = parameters.replace("[positive]", positivecheck)


            #parameters = parameters.replace("[LocalDatabasePath]", str(LocalDatabasePath))
            parameters = parameters.replace("[LocalDatabasePath]", str(filenamecandidates))


            print(dataset, accession, formula, ionform)
            parameterfile = str(softwareinput / "parameter.txt")

            with open(parameterfile, 'w', encoding="utf-8", newline="") as f:
                f.write(parameters)

            #
            #
            # Exec MSFINDER
            #
            print("lauching MetFrag", len(candidates_lst), "candidates")
            ddatoolbox.writeLog(logfilename, ["MetFrag", accession, "launched", len(candidates_lst)])

            #commanttext = "C:/Program Files/Java/jdk-22/bin/java.exe -jar MetFragCommandLine-2.5.0.jar path\parameter.txt"
            #commanttext = '"C:/Program Files/Java/jdk-22/bin/java.exe" -jar ./Software/MetFragCommandLine-2.5.0.jar path\parameter.txt'
            commanttext = '"{java:s}" -jar {metflag:s} {parameter:s}'
            #java = "C:/Program Files/Java/jdk-22/bin/java.exe"
            #metflag = str(softwareroot / 'MetFragCommandLine-2.6.1.jar')

            java = paths.getsettings("SOFTWAREPATH", "java")
            execmetfrag = paths.getsettings("SOFTWAREPATH", "execmetfrag")

            commanttext = commanttext.format(java = java, metflag = execmetfrag, parameter = parameterfile)

            subprocess.run(commanttext, shell=True)
            counter += 1


    scoredict = {}
    outputtable = []
    pathlist = list(softwareoutput.glob(r'*'))

    #print(paths, softwareoutput)


    for path in pathlist:
        accession = re.findall(r'(MSBNK-[A-Za-z0-9_]+-[A-Za-z0-9]+)', str(path))
        #
        # Get MassBank record ID
        #
        if len(accession) == 0:
            continue
        accession = accession[0]


        accessionfull = re.findall(r'(MSBNK-[A-Za-z0-9_]+-[A-Za-z0-9]+_[0-9]+)', str(path))
        #
        # Get MassBank record ID
        #

        accessionfull = accessionfull[0]
        ionform = ionformdict[accessionfull]
        #
        # Read the MetFrag result CSV
        #
        rows = ddatoolbox.readsmalltable(path,["Score", "InChIKey","MolecularFormula", "Name", "CompoundName"], filetype = "csv")
        # If the result is empty
        if len(rows) == 0:
            continue
        rank = 1

        for row in rows:

            score=float(row["Score"])
            rank = int(rank)
            inchikey = row["InChIKey"]
            formula=row["MolecularFormula"]
            if "CompoundName" in row:
                name = row["CompoundName"]
            else:
                name = row["Name"]

            if rank < rankingthreshold:
                if not accession in scoredict:
                    scoredict[accession] = []
                temp = {"score":score, "ionform": ionform, "name": name, "formula": formula, "InChIKey": inchikey}
                scoredict[accession].append(temp)
            rank += 1




    for accession, hits in scoredict.items():
        counter = 0
        hits_sorted = sorted(hits, key = lambda x: x["score"], reverse=True)
        for hit in hits_sorted:
            id = accession
            type = "compound"
            scoreing = "MetFrag"
            db = useddb# = "MS-FINDER_DB"

            score = hit["score"]
            inchi = hit["InChIKey"]
            name = hit["name"]
            ionform = hit["ionform"]
            outputtable.append([id, type, scoreing, db, score, inchi, ionform, name])
            counter += 1
            if counter > rankingthreshold:
                break
            #ddatoolbox.writeLog(logfilename, ["MetFrag", id, type, scoreing, db, score, inchi, ionform, name])


    with open(outputpath/ "output.txt", mode="a",encoding='UTF-8', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(outputtable)

    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Finished", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), softwarename])
    ddatoolbox.writeLog(startlogfilename, ["Finished", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), softwarename])



def specpred(parameter, containername):
    """Predict a compound's MS/MS spectra with CFM-ID and write them as MassBank records.

    Runs ``cfm-predict`` for the given compound (via a dockerized CFM-ID 4.x
    container), retrying up to 5 times on failure, and parses the resulting
    per-collision-energy (10/20/40 V) fragment spectra and SMILES fragment
    annotations. For each collision energy with a non-empty predicted spectrum,
    builds and writes out a synthetic MassBank record (accession prefixed with
    ``MSBNK-CFMIDP_``/``MSBNK-CFMIDN_``) under ``./CFM-ID/pos/`` or
    ``./CFM-ID/neg/<InChIKey block 1>/``. If CFM-ID reports the molecule as
    unsupported/unparseable, writes an ``unacceptable.txt`` marker file instead
    and stops retrying.

    Args:
        parameter (list): A 7-element list, in order: ``[smiles, mw, name,
            formula, inchikey, softwareoutput, polarity]``.

            - smiles (str): SMILES string of the compound.
            - mw (float): Monoisotopic mass of the neutral compound.
            - name (str): Compound display name (non-ASCII characters are
              stripped before use).
            - formula (str): Molecular formula.
            - inchikey (str): Full InChIKey of the compound.
            - softwareoutput: Unused directly in this function (output path is
              derived from ``predictedpath``/``inchikey``/``polarity`` instead).
            - polarity (str): ``"POSITIVE"`` or ``"NEGATIVE"``.
        containername (str): Name of the running Docker container hosting CFM-ID,
            used to target ``docker exec``/``docker cp`` commands.

    Returns:
        tuple: A tuple of five elements.

            - inchikey1st (str): First (14-character) block of the InChIKey.
            - status (str): Final status string (e.g. ``"Successfully finished"``,
              a CFM-ID error message, or ``"not processed"``).
            - containername (str): The Docker container name used (echoed back).
            - smiles (str): The input SMILES (echoed back).
            - Failcounter (int): Number of prediction attempts made (1-5).
    """
    [smiles, mw, name, formula, inchikey, softwareoutput, polarity] = parameter


    def generateunacceptable(path, status):
        with open(path, mode="w") as f:
            f.write(status)


    # Replace any non-ASCII characters with an empty string
    name = re.sub(r'[^\x00-\x7F]+', '', name)
    #print(cleaned) # output: Python2:


    inchikey1st = inchikey.split("-")[0]

    predictedpath = "./CFM-ID/"
    predictedpath = Path(predictedpath).resolve()

    if polarity == "POSITIVE":
        ionform = "[M+H]+"
        precmz = float(mw) + 1.0072764
        outputpath = predictedpath / "pos" / inchikey1st

        accessionf = "MSBNK-CFMIDP_{:s}_{:d}V"

    else:
        ionform = "[M-H]-"
        precmz = float(mw) - 1.0072764
        outputpath = predictedpath / "neg" / inchikey1st
        accessionf = "MSBNK-CFMIDN_{:s}_{:d}V"


    Failcounter = 0
    status = "not processed"

    if not os.path.isdir(outputpath):
        os.mkdir(outputpath)


    while(Failcounter < 5):

        Failcounter += 1

        command_prediction = 'docker exec -it {cname:s} sh -c "cd /cfmid/public/; cfm-predict \'{smiles:s}\' 0.001 /trained_models_cfmid4.0/{ionform:s}/param_output.log /trained_models_cfmid4.0/{ionform:s}/param_config.txt 1 output.txt"'
        command_prediction = command_prediction.format(ionform = ionform, cname= containername, smiles = smiles)
        result = subprocess.run(command_prediction, shell=True, capture_output=True, text=True)

        time_duration = 0.5
        time.sleep(time_duration)


        command_output_data_copy = 'docker cp {cname:s}:"/cfmid/public/output.txt" "{output:s}"'.format(output = str(outputpath), cname= containername)
        resultcopy = subprocess.run(command_output_data_copy, shell=True, capture_output=True, text=True)

        time.sleep(time_duration)

        spec = {}
        spec['energy0'] = {}
        spec['energy1'] = {}
        spec['energy2'] = {}
        anno = {}
        energy = ""

        outputfilename = str(outputpath/"output.txt")
        try:
            with open(outputfilename) as fo:
                reader = csv.reader(fo, delimiter=' ')
                for row in reader:
                    #print(row)

                    if len(row) == 0:
                        energy =  "anno"
                        continue

                    if "#" in row[0]:
                        continue

                    if "energy" in row[0]:
                        energy =  row[0]
                        continue

                    if len(row) < 3:
                        continue

                    if energy == "anno":
                        frgid = row[0]
                        frgsmiles = row[2]
                        anno[frgid] = frgsmiles
                        continue

                    frgmz = float(row[0])
                    frgrelint = float(row[1])
                    frgid = row[2]
                    spec[energy][frgmz] = {"int": frgrelint, "frgid": frgid}
        except:
            continue


        for targetcollision in [10, 20, 40]:

            if targetcollision == 10:
                energy = "energy0"
            elif targetcollision == 20:
                energy = "energy1"
            else:
                energy = "energy2"

            if len(spec[energy]) == 0:
                continue

            specdata = {}
            annodata = {}
            for frgmz, item in spec[energy].items():
                if item["frgid"] in anno:
                    annotationsmiles = anno[item["frgid"]]
                else:
                    annotationsmiles = "nd"
                specdata[frgmz] = item["int"]
                annodata[frgmz] = {
                "annotation": annotationsmiles,
                "intensity": item["int"]
                }

            mb = ddatoolbox.massbankio.MassBankRecord()

            mb.add("AUTHORS:", "Shin-MassBank project")
            mb.add("LICENSE:", "CC BY-NC-SA")
            mb.add("PUBLICATION:", "Shin-MassBank project in preparation")
            mb.add("CH$NAME:", name)
            mb.add("CH$FORMULA:", formula)
            mb.add("CH$EXACT_MASS:", mw)
            mb.add("CH$SMILES:", smiles)
            #mb.add("CH$IUPAC:", inchi)
            mb.add2("CH$LINK:", "INCHIKEY", inchikey)
            #mb.add("CH$LINK") = "CID:11353257
            mb.add("AC$INSTRUMENT:", "Spectra Predicton by CFM-ID 4.4.7")
            mb.add("AC$INSTRUMENT_TYPE:", "LC-ESI-QFT")
            mb.add2("AC$MASS_SPECTROMETRY:","MS_TYPE", "MS2")
            mb.add2("AC$MASS_SPECTROMETRY:","ION_MODE", polarity)
            mb.add2("AC$MASS_SPECTROMETRY:","COLLISION_ENERGY", "{:s} V".format(str(targetcollision)))
            mb.add2("AC$MASS_SPECTROMETRY:","FRAGMENTATION_MODE", "CID")
            mb.add2("AC$MASS_SPECTROMETRY:","IONIZATION", "ESI")
            mb.add2("MS$FOCUSED_ION:","ION_TYPE", ionform)
            mb.add2("MS$FOCUSED_ION:","PRECURSOR_M/Z", precmz)


            mb.add_spectra_for_MassBank(specdata)
            mb.add_annotation_from_items(annodata, ["int", "annotation"])
            mb.add_recordtitle()

            accession = accessionf.format(inchikey1st, targetcollision)
            mb.add("ACCESSION:", accession)
            mb.add_date()
            mb.show(path = str(outputpath) + "/")


        filenametemp = accessionf.format(inchikey1st, 40)
        filenametemp = filenametemp+".txt"
        pathforcheck = outputpath / filenametemp
        filenametemp2 = "unacceptable.txt"
        pathforunacceptable = outputpath / filenametemp2


        status = result.stdout
        status = re.split(r'\n', status)[0]

        if os.path.isfile(pathforcheck):
            #print(filenametemp, "Successfully finished")
            status = "Successfully finished"
            break
        elif "Unsupported input molecule:" in status:
            generateunacceptable(pathforunacceptable, status)
            break
        elif "Could not ionize" in status:
            generateunacceptable(pathforunacceptable, status)
            break
        elif "Explicit valence" in status:
            generateunacceptable(pathforunacceptable, status)
            break
        elif "kekulize mol" in status:
            generateunacceptable(pathforunacceptable, status)
            break



    return(inchikey1st, status, containername, smiles, Failcounter)




def getfilelist(inputpath, accession_stem_ms2):
    """List MassBank MS2 record files matching a glob pattern and accession stem.

    Args:
        inputpath (str or Path): Glob pattern for candidate record files.
        accession_stem_ms2 (str): MS2 accession-stem substring that must appear
            in the filename (e.g. ``"HSA001P"``).

    Returns:
        list[str]: Paths to matching ``.txt`` record files.
    """
    files = glob.glob(str(inputpath))
    filelist = []
    for file in files:
        if not ".txt" in file:
            continue
        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)
    return(filelist)

def loadcompdata(softwareroot):
    """Load the MS-FINDER structure database plus its user difference-library extension.

    Reads MS-FINDER's ``MsfinderStructureDB-VS15.esd`` reference table and the
    supplementary ``difflibrary_top2.txt`` user database, each sorted by exact
    mass, for use as CFM-ID prediction candidates.

    Args:
        softwareroot (Path): Root path of the ``software`` folder, containing
            ``MSFINDER/Resources/MsfinderStructureDB-VS15.esd`` and
            ``difflibrary_top2.txt``.

    Returns:
        list: A 2-element list ``[csv_db, csv_db_diff]``.

            - csv_db (list[dict]): Rows of the MS-FINDER structure DB, each with
              keys ``"Title"``, ``"Exact mass"``, ``"Formula"``, ``"SMILES"``,
              ``"InChIkey"``, sorted by exact mass.
            - csv_db_diff (list[dict]): Rows of the user difference library, same
              keys (``"InChIkey"`` copied from the file's ``"InChIKey"`` column),
              sorted by exact mass.
    """
    #
    # Load compound db
    #
    filenametemp = softwareroot/ "MSFINDER"/"Resources"/"MsfinderStructureDB-VS15.esd"

    print(filenametemp)
    csv_db = ddatoolbox.readsmalltable(filenametemp, keywords = ["Title", "Exact mass", "Formula", "SMILES", "InChIkey"], filetype = "tab")
    csv_db = sorted(csv_db, key = lambda x: x["Exact mass"])

    #
    # added 20251207
    #
    filenametemp = softwareroot/"difflibrary_top2.txt"
    print(filenametemp)

    csv_db_diff = ddatoolbox.readsmalltable(filenametemp, keywords = ["Title", "Exact mass", "Formula", "SMILES", "InChIKey"], filetype = "tab")
    csv_db_diff = sorted(csv_db_diff, key = lambda x: x["Exact mass"])
    for csv_db_diff_temp in csv_db_diff:
        csv_db_diff_temp["InChIkey"] = csv_db_diff_temp["InChIKey"]
    return([csv_db, csv_db_diff])




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


def getcandicates(inputpath, accession_stem_ms2, polarity, ms1thres, production_correlation_alpha, csv_dbs):
    """Collect candidate compounds within precursor-mass tolerance across a dataset.

    For every record matching ``inputpath``/``accession_stem_ms2`` (skipping
    records whose correlation-cleaned spectrum has fewer than 3 peaks), computes
    the neutral compound mass implied by a single ``[M+H]+``/``[M-H]-`` adduct
    (the only forms the Docker CFM-ID image supports) and collects every entry
    from ``csv_dbs`` whose exact mass falls within the m/z tolerance, deduplicated
    by short InChIKey across the whole dataset.

    Args:
        inputpath (str or Path): Glob pattern for candidate record files.
        accession_stem_ms2 (str): MS2 accession-stem substring to filter files by.
        polarity (str): ``"POSITIVE"`` or ``"NEGATIVE"``.
        ms1thres (float): MS1 ppm matching threshold.
        production_correlation_alpha (float): Significance level passed to
            ``MassBankRecord.wash_spectra`` for spectral cleanup.
        csv_dbs (list[list[dict]]): One or more compound-database tables (as
            returned by :func:`loadcompdata`), each row with keys ``"Title"``,
            ``"Exact mass"``, ``"SMILES"``, ``"InChIkey"``, ``"Formula"``.

    Returns:
        dict: Mapping of short InChIKey (str) to a dict with keys ``"Title"``,
            ``"Exact mass"``, ``"SMILES"``, ``"InChIKey"``, ``"formula"`` — the
            union of all matching candidates found across the dataset.
    """
    filelist = getfilelist(inputpath, accession_stem_ms2)
    candidates_list = {}
    for ms2filename in filelist:
        #
        # Loading MassBank recored MS2
        #
        mb2 = ddatoolbox.massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)
        #
        # Getting information
        #
        accession = mb2.get_accession()

        #if not "0205016" in accession:
        #    continue
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
        #
        # Docker version of CFM-ID only accepct M+H or M-h
        #
        if polarity == "NEGATIVE":
            ionform = "[M-H]-"
            compoundmw = precmz + 1.0072764
            polsymbol = "-"
        else:
            ionform = "[M+H]+"
            compoundmw = precmz - 1.0072764
            polsymbol = "+"

        absolutethres = compoundmw * 0.000001 * ms1thres * 2
        #absolutethres = compoundmw * 0.000001 * ms1thres * 2
        compoundmwub = compoundmw + absolutethres
        compoundmwlb = compoundmw - absolutethres
        #
        # Spectra cleaning up by correlation coefficient to precursor
        #
        cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)
        if len(cleanedspectra) < 3:
            continue
        #
        # Preparation of candidate list
        #
        for compdb in csv_dbs:
            for datatemp in compdb:
                em = float(datatemp["Exact mass"])

                if compoundmwlb < em < compoundmwub:
                    inchikey = datatemp["InChIkey"]
                    inchikeyshort = inchikey.split("-")[0]

                    if not inchikeyshort in candidates_list:
                        candidates_list[inchikeyshort] = {
                        "Title": datatemp["Title"],
                        "Exact mass": datatemp["Exact mass"],
                        "SMILES": datatemp["SMILES"],
                        "InChIKey": datatemp["InChIkey"],
                        "formula": datatemp["Formula"],
                        }
                        #print(accession, "added", inchikeyshort, datatemp["Formula"])
                    else:
                        pass
                        #print(accession, "not added", inchikeyshort)
        print(accession, len(candidates_list))
    return(candidates_list)






def search_cfmid(dataset, path, containerlist, ncpu = 10):
    """Predict MS/MS spectra for candidate compounds with CFM-ID and search them like MassBank.

    Loads the MS-FINDER structure DB (via :func:`loadcompdata`), collects
    precursor-mass-matching candidates across the dataset (via
    :func:`getcandicates`, cached to ``candidates_list.pkl``), then predicts each
    not-yet-predicted candidate's spectrum in parallel across the given Docker
    containers (via :func:`specpred`, dispatched with ``joblib.Parallel``).
    Predicted spectra already present under the CFM-ID prediction cache
    (``./res/CFM-ID/``) or previously marked unacceptable are skipped. Finally,
    all predicted MassBank-format records for this dataset's candidates are
    copied into a local spectral-library folder and indexed for
    spectral-similarity search (via ``ddatoolbox.massbankserver``), and every
    step2 record is searched against it, with matches appended to the dataset's
    ``output.txt``.

    Args:
        dataset (str): Dataset identifier string.
        path (dict): Dictionary containing file path information. Must include
            ``"logfilename"`` and ``"startlogfilename"``; other paths are
            re-derived from ``dataset`` rather than read from ``path``.
        containerlist (list[str]): Docker container names to parallelize CFM-ID
            prediction across; determines the batch size for each
            ``joblib.Parallel`` dispatch.
        ncpu (int, optional): Number of parallel worker processes for
            ``joblib.Parallel``. Defaults to ``10``.

    Returns:
        None
    """
    software = "CFM-ID"

    #datasetpath = path["dataset"]
    #softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    #inputpath = datasetpath / "step2" / "*"
    #inputpath1 = datasetpath / "step1"
    #outputpath = path["output"]
    #softwareinput = path["softwareinput"]
    #softwareoutput = path["softwareoutput"]


    datasetpath = "../data/"+dataset+"/"
    datasetpath = Path(datasetpath).resolve()
    softwarepath = "../software/"
    softwarepath = Path(softwarepath).resolve()


    # annotation output folders, per dataset
    annotationpath = datasetpath /"annotation"
    outputpath = annotationpath / software

    softwareinput = outputpath / "file"
    softwareoutput = outputpath / "result"
    massbankdata = outputpath / "MassBank-data"
    spectra = massbankdata  / "CFMID"
    massbank = outputpath / "MassBank"

    softwareroot = softwarepath

    predictedpath = "./res/CFM-ID/"
    predictedpath = Path(predictedpath).resolve()

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"

    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")
    accession_stem = settings["accession_stem"]
    accession_stem_ms2 = dataset+"P"
    accession_stem_ms1 = dataset+"S"
    polarity = settings['polarity']
    ms1thres = float(settings['ms1_mz_averaged_thres'])
    ms2thres = float(settings['ms2_mz_averaged_thres'])


    if polarity == "POSITIVE":
        ionform = "[M+H]+"
    else:
        ionform = "[M-H]-"


    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "CFM-ID"], mode = "a")
    ddatoolbox.writeLog(startlogfilename, ["Started", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "CFM-ID"], mode = "a")

    print(dataset, "loading compound data")

    csv_dbs = loadcompdata(softwareroot)
    print(dataset, "Gathering Compound InChiKey for searching")

    filename = str(outputpath / 'candidates_list.pkl')
    if not os.path.isfile(filename):
        #
        # Load compound db
        #
        candidates_list = getcandicates(inputpath, accession_stem_ms2, polarity, ms1thres, production_correlation_alpha, csv_dbs)
        print(dataset, len(candidates_list))
        with open(filename, 'wb') as f:
            pickle.dump(candidates_list, f)


    print(dataset, "Predicting spectra by CFM-ID")
    with open(filename, 'rb') as f:
        candidates_list = pickle.load(f)

    totalcandidatenum = len(candidates_list)
    candicatecounter = 0

    joblist = []
    #compoundlibrarypath = "../software/MsfinderStructureDB-VS15_forMetFrag.csv"
    for csv_db in csv_dbs:
        for datatemp in csv_db:
            #inchi = row[1]
            #mw = row[2]
            #formula = row[3]
            #inchikey1st = row[4]
            #inchikey2nd = row[5]
            #smiles = row[6]
            #name = row[7]
            #inchikey3rd = row[8]
            #inchikey = "{:s}-{:s}-{:s}".format(inchikey1st, inchikey2nd, inchikey3rd)
            inchikey = datatemp["InChIkey"]

            mw = datatemp["Exact mass"]
            formula = datatemp["Exact mass"]
            inchikey1st = inchikey.split("-")[0]
            smiles = datatemp["SMILES"]
            name = datatemp["Title"]

            #inchikey = "{:s}-{:s}-{:s}".format(inchikey1st, inchikey2nd, inchikey3rd)


            if smiles == "SMILES":
                continue
            if not inchikey1st in candidates_list:
                continue
            if polarity == "POSITIVE":
                specoutputpath = predictedpath / "pos" / inchikey1st
                accessionf = "MSBNK-CFMIDP_{:s}_40V.txt"

            else:
                specoutputpath = predictedpath / "neg" / inchikey1st
                accessionf = "MSBNK-CFMIDN_{:s}_40V.txt"


            candicatecounter += 1
            #if candicatecounter > 10:
            #    break



            filenametemp = accessionf.format(inchikey1st)
            pathforcheck = specoutputpath / filenametemp
            #print(pathforcheck)

            filenametemp = "unacceptable.txt"
            pathforunacceptable = specoutputpath / filenametemp
            if os.path.isfile(pathforcheck):
                #print(inchikey1st, "Already it exists in CFM-ID")
                #print(dataset, candicatecounter,"/",totalcandidatenum, "*", inchikey1st, name)
                continue
            elif os.path.isfile(pathforunacceptable):
                continue
            else:
                print(dataset, candicatecounter,"/",totalcandidatenum, " ", inchikey1st, name)



            joblist.append([smiles, mw, name, formula, inchikey, softwareoutput, polarity])

            if len(joblist) < len(containerlist):
                continue
            #result = joblib.Parallel(n_jobs=ncpu)(joblib.delayed(f)(optimizer.x[:,i], parameters, i) for i in range(population_size))


            result = joblib.Parallel(n_jobs=ncpu)(joblib.delayed(specpred)(parameter, containername) for parameter, containername in zip(joblist, containerlist))
            for (inchikey, state, containername, smiles, Failcounter) in result:
            #for parameter, containername in zip(joblist, containerlist):
                #inchikey, state = specpred(parameter, containername)
                print(inchikey, Failcounter, state, containername, smiles)
                if not "Successfully" in state:
                    ddatoolbox.writeLog(logfilename, ["Failed", inchikey, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), state, containername, smiles], mode = "a")
            joblist = []

    result = joblib.Parallel(n_jobs=ncpu)(joblib.delayed(specpred)(parameter, containername) for parameter, containername in zip(joblist, containerlist))
    for (inchikey, state, containername, smiles, Failcounter) in result:
    #for parameter, containername in zip(joblist, containerlist):
        #inchikey, state = specpred(parameter, containername)
        print(inchikey, Failcounter, state, containername, smiles)
        if not "Successfully" in state:
            ddatoolbox.writeLog(logfilename, ["Failed", inchikey, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), state, containername, smiles], mode = "a")
    #
    # "File Copy and creating"
    #
    print(dataset, "File copy")

    if not os.path.isdir(massbankdata):
        os.mkdir(massbankdata)
        os.mkdir(spectra)
        os.mkdir(massbank)

        totalinchikey = len(candidates_list)
        inchikeyconter = 0

        for inchikey1st in candidates_list:

            progress_bar(inchikeyconter, totalinchikey, bar_length=40)


            inchikeyconter += 1


            if polarity == "POSITIVE":
                outputpathtemp = predictedpath / "pos" / inchikey1st / "*"
            else:
                outputpathtemp = predictedpath / "neg" / inchikey1st / "*"

            files = glob.glob(str(outputpathtemp))
            filelist = []
            for file in files:
                if not ".txt" in file:
                    continue
                if not "MSBNK" in file:
                    continue
                filelist.append(file)
            for ms2filename in filelist:
                shutil.copy(ms2filename, spectra)


        print()
        print(dataset, "Generating Spec Search library")
        path = str(massbankdata/"*")
        targetpath = str(massbank) + "/"
        massbankserver.readMassBankfiles(path, targetpath, digit = 1, ionization = "ESI", polarity = polarity, MSTYPE = "MS2")

        shutil.rmtree(spectra)

    print()
    print(dataset, "Searching the spectral library predicted by CFM-ID")


    with open(str(massbank/"MassBankRecords.pkl"), "rb") as tf:
        new_dict = pickle.load(tf)

    outputtextpath = outputpath/"output.txt"

    #
    # Remove old output.txt file
    #
    if os.path.isfile(outputtextpath):
        os.remove(outputtextpath)

    filelist = getfilelist(inputpath, accession_stem_ms2)
    totalfile = len(filelist)
    fileconter = 0
    for ms2filename in filelist:


        progress_bar(fileconter, totalfile, bar_length=40)
        fileconter += 1
        #
        # Loading MassBank recored MS2
        #
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)
        #
        # Getting information
        #
        accession = mb2.get_accession()

        #if not "0205016" in accession:
        #    continue
        cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)
        if len(cleanedspectra) < 3:
            continue

        query = cleanedspectra
        precmz, stdev = mb2.get_precmz()

        absolutethres = precmz * 0.000001 * ms1thres * 2
        precmzlb = precmz - absolutethres
        precmzub = precmz + absolutethres

        output = massbankserver.searchMassBank(new_dict, query, scorethreshold = 0.0, mzthreshold = ms2thres, digit = 1, path = str(massbank) +"/", unit = "ppm")
        resultsummary = {}
        for mbaccession, datatemp in output.items():
            precmztemp = datatemp['PRECURSOR_M/Z']
            if not precmzlb < precmztemp < precmzub:
                continue
            inchikey = datatemp['INCHIKEY']
            inchikey1st = inchikey.split("-")[0]
            score = datatemp['dice']
            name = datatemp['RECORD_TITLE']
            formula = datatemp['FORMULA']
            if not inchikey1st in resultsummary:
                resultsummary[inchikey1st] = {
                "sumscore": 0,
                "scorelist": [],
                "name": name,
                "formula": formula
                }
            resultsummary[inchikey1st]["sumscore"] += score
            resultsummary[inchikey1st]["scorelist"].append(score)

        inchikey1stlist = resultsummary.keys()
        inchikey1stlist_sorted = sorted(inchikey1stlist, key = lambda x: resultsummary[x]["sumscore"], reverse= True)

        outputtext = []
        for inchikey1st in inchikey1stlist_sorted:
            temp = resultsummary[inchikey1st]
            compname = temp["name"].split("; ")[0]
            outputtemp = [accession, "compound", "CFM-ID", "MS-FINDER_DB", temp["sumscore"], inchikey1st, ionform, compname]
            outputtext.append(outputtemp)

        with open(outputpath/"output.txt", mode="a",encoding='UTF-8', newline="") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(outputtext)



    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Finished", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "CFM-ID"])
    ddatoolbox.writeLog(startlogfilename, ["Finished", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), "CFM-ID"])



def search_similarity(dataset, path, database = "GNPS"):
    """Search each record's cleaned MS2 spectrum against a GNPS or MassBank reference library.

    Loads the binary spectral-search database (via
    ``ddatoolbox.massbankserver.setupMassBank``) for the requested reference
    library and the record's polarity, then for every step2 record with at least
    one candidate formula (from ``COMMENT: PROPOSED_MF_FROM_PRECMZ``), searches
    its correlation-cleaned spectrum (further filtered to intensity >= 9 when it
    has more than 30 peaks) against the library via
    ``ddatoolbox.massbankserver.searchMassBank``. Hits whose formula matches one of
    the record's candidate formulas are aggregated per short InChIKey (summing
    hit counts, tracking per-hit cosine score/Dice/entropy-similarity), and the
    best-scoring hits are appended to the dataset's ``output.txt`` as both
    ``SIMILARITY`` (cosine-product) and ``Entropy`` (entropy-similarity) rows.

    Args:
        dataset (str): Dataset identifier, used to build the accession-stem
            prefix for locating step2 record files.
        path (dict): Dictionary containing file path information with keys
            ``"dataset"``, ``"software"``, ``"logfilename"``,
            ``"startlogfilename"``, ``"output"``, ``"softwareinput"``, and
            ``"softwareoutput"``.
        database (str, optional): Reference spectral library to search,
            ``"GNPS"`` or ``"MassBank"``. Defaults to ``"GNPS"``.

    Returns:
        None
    """
    processmassbankdata = False
    #
    # database = "msfinder"

    datasetpath = path["dataset"]
    softwareroot = path["software"]
    logfilename = path["logfilename"]
    startlogfilename = path["startlogfilename"]

    inputpath = datasetpath / "step2" / "*"
    inputpath1 = datasetpath / "step1"
    outputpath = path["output"]
    softwareinput = path["softwareinput"]
    softwareoutput = path["softwareoutput"]
    #
    accession_stem_ms2 = dataset+"P"
    accession_stem_ms1 = dataset+"S"

    if database == "GNPS":
        softwarename = "SIMILARITY"
        database = "GNPS"
    elif database == "MassBank":
        softwarename = "SIMILARITY"
        database = "MassBank"

    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Started", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), softwarename], mode = "a")
    ddatoolbox.writeLog(startlogfilename, ["Started", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), softwarename], mode = "a")
    #
    #
    #


    settings = massbankio.read_RecordTemplate(datasetpath/"RecordTemplate.txt")
    polarity = settings['polarity']
    #LocalDatabasePath = softwareroot/"kegg_2017-07-23.csv"
    if database == "GNPS":
        if polarity == "NEGATIVE":
            path = "./res/spectra/GNPSnegative/"
        else:
            path = "./res/spectra/GNPSpositive/"
    elif database == "MassBank":
        if polarity == "NEGATIVE":
            path = "./res/spectra/MassBanknegative/"
        else:
            path = "./res/spectra/MassBankpositive/"

    ionformdata = massbankio.read_ionformlist(filename = "ion_forms.txt", polarity = polarity)

    mbdict = massbankserver.setupMassBank(path)

    outputtable = []

    filelist = []

    files = glob.glob(str(inputpath))

    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        filelist.append(file)
    totalnumber = len(filelist)
    counter = 0

    for ms2filename in filelist:
        counter += 1
        #if not "PCE000173" in ms2filename:
        #    continue
        #
        # Loading MassBank recored MS2
        #
        mb2 = ddatoolbox.massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)
        #
        # Getting information
        #
        accession = mb2.get_accession()
        #precmz, intensity = mb2.get_precmz()
        try:
            precmz, intensity = mb2.get_precmz()
        except:
            print("no prec data", accession)
            continue

        ms1thres = float(mb2.get_anything("COMMENT:","MS1_M/Z_AVERAGED_THRES"))
        ms2thres = float(mb2.get_anything("COMMENT:","MS2_M/Z_AVERAGED_THRES"))

        polarity = mb2.get_polarity()

        if polarity == "NEGATIVE":
            polsymbol = "-"
            positivecheck = "False"
        else:
            polsymbol = "+"
            positivecheck = "True"

        if processmassbankdata == True:
            #
            #PROPOSED_ION_FORM
            #
            Ion_forn_list_temp = mb2.get_anything("COMMENT:", "PROPOSED_MF_FROM_PRECMZ")
            Ion_forn_list = {}
            ionform_list = {}

            for Ion_forn_temp in Ion_forn_list_temp:
                mf = Ion_forn_temp.split(" ")[0]
                ionformtemp = Ion_forn_temp.split(" ")[2]
                ionform = "["+ionformtemp+"]"+polsymbol
                ionform = ionform.replace(";", "")
                Ion_forn_list[mf] = ionform
                ionform_list[ionformtemp] = 1

            if len(Ion_forn_list) == 0:
                print("No mf", accession)
                continue

            #
            # Molecular formula searching of precursor ion
            #
            threshold_temp = ms1thres
            hitmfs = []
            for ion_form in ionformdata.keys():
                diff_temp  = ionformdata[ion_form]["mass"]
                diff_mult  = ionformdata[ion_form]['mult']

                if ion_form in ionform_list:
                    mz_Na = (precmz - diff_temp)/diff_mult
                    result, result_list = mfsearch.mfsearch(mz_Na, threshold = threshold_temp, ihdtype = "integer", Climit = 99, Hlimit = 999, Nlimit = 5, Olimit = 99, Plimit = 2, Slimit = 2, polarity = "None")
                    difflist = ["{:6f}".format(mfsearch.calcmw(mf) - mz_Na)  for mf in result_list]

                    hitlist = []
                    for (formula_temp, diff_temp) in zip(result_list, difflist):
                        hitlist.append((formula_temp, diff_temp))

                    for hitmf, hitdiff in zip(result_list, difflist):
                        hitdiffppm = float(hitdiff)/mz_Na * 1000000

                        if polarity == "POSITIVE":
                            ionform = "[{:s}]+".format(ion_form)
                        else:
                            ionform = "[{:s}]-".format(ion_form)

                        texttemp = "{} of {} with delta m/z {:4f} ppm".format(mfsearch.cleanmf(hitmf), ion_form, hitdiffppm)
                        mb2.add2("COMMENT:","PROPOSED_MF_FROM_PRECMZ", texttemp)

        #
        #PROPOSED_ION_FORM
        #
        Ion_forn_list_temp = mb2.get_anything("COMMENT:", "PROPOSED_MF_FROM_PRECMZ")

        Ion_forn_list = {}
        for Ion_forn_temp in Ion_forn_list_temp:
            mf = Ion_forn_temp.split(" ")[0]
            ionform = Ion_forn_temp.split(" ")[2]
            ionform = "["+ionform+"]"+polsymbol
            ionform = ionform.replace(";", "")
            Ion_forn_list[mf] = ionform

        if len(Ion_forn_list) == 0:
            print("No mf", accession)
            continue

        #
        # Spectra cleaning up
        #
        cleanedspectra_temp = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)
        if len(cleanedspectra_temp) < 3:
            continue

        cleanedspectra = {}

        for mz, intensity in cleanedspectra_temp.items():
            if len(cleanedspectra_temp)> 30 and intensity < 9:
                continue
            cleanedspectra[mz] = intensity



        results = massbankserver.searchMassBank(mbdict, cleanedspectra, scorethreshold = 0.01, mzthreshold = 0.02, digit = 2, path = path, unit = "Da")

        hittedMassBankrecordid = {}

        for sid, hittedspectra in results.items():
            #if not sid == "CCMSLIB00005463808":
            #    continue
            hittedformula = mfsearch.cleanmf(hittedspectra["FORMULA"])
            hittedthormz = mfsearch.calcmw(hittedformula, polarity)
            recordtitle = hittedspectra["RECORD_TITLE"]
            score = hittedspectra["score"]
            dice = hittedspectra["dice"]
            inchikey = hittedspectra["INCHIKEY"]
            if inchikey == "":
                continue

            if processmassbankdata == True:
                if float(score) > 0.99:
                    continue

            inchikey_1st = inchikey.split("-")[0]
            #print(hittedformula, Ion_forn_list)
            if hittedformula in Ion_forn_list:
                ionform = Ion_forn_list[hittedformula]
                #print(hittedspectra)

                spectratemp = hittedspectra["spectra"]
                specent = calcsimilarity(cleanedspectra, spectratemp, 0.02)
                #print(sid, specent)


                if not inchikey_1st in hittedMassBankrecordid:
                    hittedMassBankrecordid[inchikey_1st]={
                    "recordtitle": recordtitle,
                    "score":[score],
                    "dice":[dice],
                    "specent":[specent],
                    "sid": sid,
                    "mf":hittedformula,
                    "counter": 1,
                    "ionform": ionform,
                    "name": recordtitle
                    }
                else:
                    hittedMassBankrecordid[inchikey_1st]["counter"] += 1
                    hittedMassBankrecordid[inchikey_1st]["score"].append(score)
                    hittedMassBankrecordid[inchikey_1st]["dice"].append(dice)
                    hittedMassBankrecordid[inchikey_1st]["specent"].append(specent)
        #
        # Sort result
        #

        #for inchikey in sorted(hittedMassBankrecordid.keys(), key=lambda x: hittedMassBankrecordid[x]["counter"] + hittedMassBankrecordid[x]["score"], reverse=True):
        for inchikey in sorted(hittedMassBankrecordid.keys(), key=lambda x: max(hittedMassBankrecordid[x]["score"]), reverse=True):

            recordtitle=hittedMassBankrecordid[inchikey]["recordtitle"]
            score=hittedMassBankrecordid[inchikey]["score"]
            dice=hittedMassBankrecordid[inchikey]["dice"]
            sid=hittedMassBankrecordid[inchikey]["sid"]
            mf=hittedMassBankrecordid[inchikey]["mf"]
            hitcounter= hittedMassBankrecordid[inchikey]["counter"]
            ionform = hittedMassBankrecordid[inchikey]["ionform"]
            name = hittedMassBankrecordid[inchikey]["name"]
            specent = hittedMassBankrecordid[inchikey]["specent"]

            #print(dataset, "MassBank", inchikey, hitcounter, score, mf, sid, recordtitle)
            #outputtemp = [accession, "compound", softwarename, database, hitcounter, inchikey, ionform,name]

            sumscore = sum(score)
            sumdice = sum(dice)

            score = max(score)
            dice = max(dice)
            specent = max(specent)

            nameplus = "{:s}::{:s}::{:d}::{:.4f}::{:.4f}::{:.4f}".format(name, mf, hitcounter, score, dice, specent)
            outputtemp = [accession, "compound", softwarename, database, score, inchikey, ionform, nameplus]
            print(counter, "/", totalnumber, outputtemp)
            outputtable.append(outputtemp)


            outputtemp = [accession, "compound", "Entropy", database, specent, inchikey, ionform, nameplus]
            print(counter, "/", totalnumber, outputtemp)
            outputtable.append(outputtemp)


    with open(outputpath/ "output.txt", mode="a",encoding='UTF-8', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(outputtable)

    dt_now = datetime.datetime.now()
    ddatoolbox.writeLog(logfilename, ["Finished", accession_stem_ms2, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), softwarename])
    ddatoolbox.writeLog(startlogfilename, ["Finished", "annotation", dt_now.strftime('%Y.%m.%d.%H.%M.%S'), softwarename])

def calcentoropy(spec):
    """Compute the (weighted) Shannon entropy of a spectrum, per Li et al.'s spectral entropy method.

    Normalizes ``spec`` to unit total intensity and computes its Shannon entropy.
    When the entropy is below 3 (the low-entropy regime where entropy-based
    similarity becomes unreliable per Li et al., 2021, Nature Methods), applies
    the recommended intensity-reweighting correction
    (``intensity ** (0.25 * (1 + entropy))``) and recomputes the entropy on the
    reweighted, renormalized spectrum.

    Args:
        spec (dict): Input spectrum, mapping m/z to intensity.

    Returns:
        tuple or int: A tuple ``(entropy, weighted_spectrum)`` where ``entropy``
            (float) is the (possibly weight-corrected) Shannon entropy and
            ``weighted_spectrum`` (dict) is the corresponding normalized spectrum
            used to compute it (only reweighted when the low-entropy correction
            was applied). Returns ``0`` (a bare int) if ``spec`` has zero total
            intensity.
    """
    spec = normalizespec(spec)

    totalintensity = sum([spec[x] for x in spec.keys()])
    if totalintensity == 0:
        return(0)

    s = 0
    spectemp = {}
    for mz, intensity in spec.items():
        if intensity <= 0:
            continue
        Itemp = intensity/totalintensity
        stemp = Itemp * math.log(Itemp)
        spectemp[mz] = Itemp
        s -= stemp

    if s < 3:
        w = 0.25 * (1 + s)
        s = 0
        spectemp = {}
        for mz, intensity in spec.items():
            if intensity <= 0:
                continue
            Itemp = (intensity/totalintensity)**w
            spectemp[mz] = Itemp
        spectemp = normalizespec(spectemp)
        for mz, intensity in spectemp.items():
            stemp = intensity * math.log(intensity)
            s -= stemp


    return(s,  spectemp)


def calcentoropy_nocheck(spec):
    """Compute the plain Shannon entropy of a spectrum, without the low-entropy correction.

    Simpler variant of :func:`calcentoropy` that always computes the unweighted
    Shannon entropy, skipping the intensity-reweighting correction applied there
    for low-entropy spectra.

    Args:
        spec (dict): Input spectrum, mapping m/z to intensity.

    Returns:
        tuple or int: A tuple ``(entropy, normalized_spectrum)``, or ``0`` (a
            bare int) if ``spec`` has zero total intensity.
    """
    spec = normalizespec(spec)

    totalintensity = sum([spec[x] for x in spec.keys()])
    if totalintensity == 0:
        return(0)
    s = 0
    for mz, intensity in spec.items():
        if intensity <= 0:
            continue
        Itemp = intensity/totalintensity
        stemp = Itemp * math.log(Itemp)
        s -= stemp

    return(s,  normalizespec(spec))


def normalizespec(spec):
    """Normalizes the intensities of a spectrum so that the sum equals 1.0.

    Args:
        spec (dict): Input spectrum {m/z: intensity}.

    Returns:
        dict: Normalized spectrum.
    """
    sumintensity = 0
    for mz, intensity in spec.items():
        sumintensity += intensity
    outputspec = {}

    for mz, intensity in spec.items():
        if intensity <= 0:
            continue
        relint = intensity/sumintensity
        outputspec[mz] = relint

    return(outputspec)


def generage_compositespectra(spec1, spec2, thres):
    """Generates a composite spectrum from two spectra by averaging matching peaks.

    Args:
        spec1 (dict): First spectrum.
        spec2 (dict): Second spectrum.
        thres (float): Tolerance for m/z matching.

    Returns:
        dict: Composite spectrum.
    """
    mz1sorted = sorted([x for x in spec1.keys()], key = lambda x: spec1[x], reverse= True)
    mz2sorted = sorted([x for x in spec2.keys()], key = lambda x: spec2[x], reverse= True)

    compsite = {}

    for mz1 in mz1sorted:
        for mz2 in mz2sorted:
            if abs(mz1-mz2) < thres:
                if not mz2 in compsite:
                    compsite[mz2] = (spec1[mz1] + spec2[mz2])/2
                    #print(mz1, mz2, compsite[mz2] )
                    break
        else:
            compsite[mz1] = spec1[mz1]/2

    for mz2 in mz2sorted:
        if not mz2 in compsite:
            compsite[mz2] = spec2[mz2]/2

    return(compsite)



def calcsimilarity(spec1, spec2, thres):
    """Compute the entropy-based spectral similarity between two spectra.

    Implements Li et al.'s (2021, Nature Methods) entropy-similarity measure:
    computes the (weight-corrected) Shannon entropy of each spectrum
    individually (via :func:`calcentoropy`) and of their merged composite
    spectrum (via :func:`generage_compositespectra` and
    :func:`calcentoropy_nocheck`), then combines them into a normalized
    similarity score in ``[0, 1]``.

    Args:
        spec1 (dict): First spectrum, mapping m/z to intensity.
        spec2 (dict): Second spectrum, mapping m/z to intensity.
        thres (float): m/z matching tolerance used to merge peaks between the
            two spectra's entropy-weighted forms.

    Returns:
        float: Entropy-similarity score, where 1.0 indicates identical spectra.
    """
    spec1 = normalizespec(spec1)
    SA, specA  = calcentoropy(spec1)
    spec2 = normalizespec(spec2)
    SB, specB  = calcentoropy(spec2)
    compspec = generage_compositespectra(specA, specB, thres)
    SAB, specAB = calcentoropy_nocheck(compspec)
    similarity = 1-(2*SAB-SA-SB)/math.log(4)
    return(similarity)

