#-------------------------------------------------------------------------------
# Name:        DDAtools
# Purpose:     Miscellaneous command-line utilities for the ShinMassBank pipeline:
#              .msp spectral-library export, pre-processing/annotation quality
#              statistics, exact-mass table generation, and conversion of the
#              MS-FINDER structure database for use with MetFrag.
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import ddatoolbox
import csv, itertools, sys, glob
import os, re,  datetime, statistics
import argparse, pickle

def generage_msp(datasets):
    """Export MassBank records as NIST-style .msp spectral-library files.

    For each dataset, reads every step4 MassBank record, determines its
    confidence/curation level (``L5``, ``L4``, ``L3-SP``, or ``L3-CL``, from the
    ``COMMENT: CONFIDENCE`` field) to decide which fields (name, formula, InChIKey,
    SMILES) are populated, then formats the record as an .msp entry (``NAME``,
    ``PRECURSORMZ``, ``FORMULA``, ``Ontology``, ``INCHIKEY``, ``SMILES``,
    ``RETENTIONTIME``, ``CCS``, ``IONMODE``, ``INSTRUMENTTYPE``, ``INSTRUMENT``,
    ``COLLISIONENERGY``, ``Comment``, ``Num Peaks``, followed by the peak list).
    Entries are grouped by ion polarity and written to one ``.msp`` file per
    polarity under the dataset's ``msp`` output path.

    Args:
        datasets (list[str]): Dataset identifiers to export.

    Returns:
        None

    Note:
        This function references a module-level ``filtered`` variable (expected to
        be ``"filtered"`` or something else) that is not defined as a parameter or
        local variable; it must be set in the caller's scope (e.g. interactively)
        before this function is called, otherwise a ``NameError`` will be raised.
    """
    targettissue = "all"

    output = {} 
    output["POSITIVE"] = []  
    output["NEGATIVE"] = []  

    for dataset in datasets:
        #dataset = "HSA{:03}".format(int(i))
        #path = f"../data/{dataset}/"
        #if not os.path.exists(path):
        #    continue
        #
        #print(dataset)

        paths = ddatoolbox.paths.generate_paths_data(dataset)

        filelist =paths["step4inputlist"]
        finalpath = paths["finalpath"]
        outputpath = paths["msp"]


        settings = ddatoolbox.massbankio.read_RecordTemplate_dataset(dataset)
        accession_stem = settings["accession_stem"]
        accession_stem_ms2 = accession_stem+"P"
        accession_stem_ms1 = accession_stem+"S"

        counter = 0


        for ms2filename in filelist:
            counter += 1
            #if counter > 10:
            #    break
            #

            mb2 = ddatoolbox.massbankio.MassBankRecord()
            mb2.read_MassBankRecord(ms2filename)

            names = mb2.get_anything("RECORD_TITLE:", "*").split("; ")
            

            confidences = mb2.get_anything("COMMENT:", "CONFIDENCE")[0].split("; ")
            confidence = confidences[0]

            if "L5" in confidence:

                level = "L5"
                #FORMULA: C39H42O2
                formula = mb2.get_formula()
                #INCHIKEY: 
                inchikey = ""
                #SMILES: CCCC#Cc1cc(C#CCCC)cc(C#Cc2cc(C#CC(C)(C)C)c(cc2C#CC(C)(C)C)C2OCCO2)c1
                smiles =  mb2.get_smiles()
                name = "Uncharacterized_"+level
                

            if "L4" in confidence:
                level = "L4"
                #FORMULA: C39H42O2
                formula = confidences[1]
                inchikey = ""
                smiles =  mb2.get_smiles()
                name = names[0]+"_"+level

            if "L2" in confidence:
                level = "L2"
                #FORMULA: C39H42O2
                formula = mb2.get_formula()
                inchikey = mb2.get_inchiley()
                smiles =  mb2.get_smiles()
                name = names[0]+"_"+level

            if "L3-SP" in confidence:
                level = "L3-SP"
                #FORMULA: C39H42O2
                formula = mb2.get_formula()
                inchikey = mb2.get_inchiley()
                smiles =  mb2.get_smiles()
                name = names[0]+"_"+level

            if "L3-CP" in confidence:
                level = "L3-CP"
                #FORMULA: C39H42O2
                if len(confidences) > 3:
                    formula = confidences[2]
                else:
                    formula = mb2.get_formula()
                inchikey = ""
                smiles =  mb2.get_smiles()
                name = names[0]+"_"+level


            #PRECURSORMZ: 81.520546
            precmz, sd = mb2.get_precmz()
            #FORMULA: C39H42O2
            #Ontology: MSBNK-SMB_Measured-HSA015P0081000
            accession = mb2.get_accession()

            #RETENTIONTIME: 13.137
            rt,sd = mb2.get_rt()
            #CCS: 
            ccs = ""
            #IONMODE: POSITIVE
            polarity = mb2.get_polarity()
            #INSTRUMENTTYPE: LC-ESI-QFT
            instrumenttype = mb2.get_instrument_type()
            #INSTRUMENT: Vanquish Flex UHPLC system, Thermo Fisher Scientific - Orbitrap Q Exactive Plus mass spectrometer, Thermo Fisher Scientific
            instrument = mb2.get_anything("AC$INSTRUMENT:", "*")
            #COLLISIONENERGY: Ramped collision energies of 40 V
            ce = mb2.get_anything("AC$MASS_SPECTROMETRY:", "COLLISION_ENERGY")

            #Comment: MSBNK-SMB_Measured-HSA015P0081000
            tissue = mb2.get_anything("SP$SAMPLE:", "*")


            
            #washedspec = mb2.wash_spectra(alpha=0.01, intthres= 0)
            #mb2.add_spectra(washedspec)
            #if filtered == "filtered":
            washedspec = mb2.wash_spectra(alpha=0.01, intthres= 0)
            mb2.add_spectra(washedspec)

            spectra = mb2.get_spectra()
            #Num Peaks: 14
            numpeaks = len(spectra)
            if numpeaks == 0:
                print("skipped", accession)



            #53.038597	77

            """
            NAME: Uncharacterized
            PRECURSORMZ: 81.520546
            FORMULA: C39H42O2
            Ontology: MSBNK-SMB_Measured-HSA015P0081000
            INCHIKEY: 
            SMILES: CCCC#Cc1cc(C#CCCC)cc(C#Cc2cc(C#CC(C)(C)C)c(cc2C#CC(C)(C)C)C2OCCO2)c1
            RETENTIONTIME: 13.137
            CCS: 
            IONMODE: POSITIVE
            INSTRUMENTTYPE: LC-ESI-QFT
            INSTRUMENT: Vanquish Flex UHPLC system, Thermo Fisher Scientific - Orbitrap Q Exactive Plus mass spectrometer, Thermo Fisher Scientific
            COLLISIONENERGY: Ramped collision energies of 40 V
            Comment: MSBNK-SMB_Measured-HSA015P0081000
            Num Peaks: 14
            53.038597	77
            54.033838	70
            55.054240	159
            56.964779	109
            64.978649	89
            67.509814	21
            79.054190	121
            80.049454	999
            81.044798	32
            81.052788	59
            81.069848	138
            82.989258	168
            83.049179	79
            83.060349	94
            """
            thisdata = []
            thisdata.append([])
            thisdata.append([f"NAME: {name}"])

            thisdata.append([f"PRECURSORMZ: {precmz}"])
            thisdata.append([f"FORMULA: {formula}"])
            thisdata.append([f"Ontology: {accession}"])
            thisdata.append([f"INCHIKEY: {inchikey}"])
            thisdata.append([f"SMILES: {smiles}"])
            thisdata.append([f"RETENTIONTIME: {rt}"])
            thisdata.append([f"CCS: {ccs}"])
            thisdata.append([f"IONMODE: {polarity}"])
            thisdata.append([f"INSTRUMENTTYPE: {instrumenttype}"])
            thisdata.append([f"INSTRUMENT: {instrument}"])
            thisdata.append([f"COLLISIONENERGY: {ce}"])
            thisdata.append([f"Comment: {level}; {tissue}"])
            thisdata.append([f"Num Peaks: {numpeaks}"])

            for mz, intensity in spectra.items():
                thision = [mz, intensity]
                thisdata.append(thision)
            

            output[polarity].extend(thisdata)

    for polarity in ["POSITIVE", "NEGATIVE"]:
        filename = f"smb_{targettissue}_{filtered}_{polarity}.msp"
        thisoutput =outputpath +"/"+filename
        with open(thisoutput, "w", encoding="UTF-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(output[polarity])        

########################################################################
def statistics_preprocessing(datasets, spec_option="", keyword=""):
    """Summarize spectral-preprocessing (averaging/filtering) statistics per dataset.

    For each dataset, reads the ``averaging`` step's log file to collect raw/used/
    averaged spectrum and peak counts, computes their per-spectrum averages, reads
    the dataset's ``<dataset>_mzcomparison.csv`` file to compute mean absolute
    precursor/product m/z deviations (raw, averaged, and calibrated-averaged), and
    re-derives the correlation-filtered peak/spectrum counts by re-running
    ``MassBankRecord.wash_spectra`` over the step4 record files. Results for every
    dataset are printed, written to ``statisticsprocessing.txt``, and appended to
    the shared MassBank log file.

    Args:
        datasets (list[str]): Dataset identifiers to process (currently unused;
            see the Note below).
        spec_option (str, optional): Reserved for future spectral-processing
            options; currently unused. Defaults to ``""``.
        keyword (str, optional): Filter string; only step4 record files whose
            filename contains ``keyword`` are processed. If empty, all files are
            processed. Defaults to ``""``.

    Returns:
        None

    Note:
        This function iterates over a module-level ``datasettemp`` variable rather
        than the ``datasets`` parameter; ``datasettemp`` must be set in the
        caller's scope before this function is called, otherwise a ``NameError``
        will be raised. It also uses ``datetime`` without importing the
        ``datetime`` module in this file, which will raise a ``NameError`` unless
        ``datetime`` has been made available some other way (e.g. via ``ddatoolbox``
        re-exporting it into the global namespace).
    """
    dt_now = datetime.datetime.now()
    logfilename = "../records/log/MassBank"+ dt_now.strftime('%Y.%m.%d') + ".txt"

    targets = [
        "spec_raw",  # total
        "spec_used",  # total
        "spec_averaged",  # total
        "spec_filtered",  # total
        "peak_raw",  # total
        "peak_used",  # total
        "peak_averaged",  # total
        "peak_filtered",  # total
        "ave_peak_raw",  # total
        "ave_peak_used",  # total
        "ave_peak_averaged",  # total
        "ave_peak_filtered",  # total

        "diff_prec_avecalib",
        "diff_prec_ave",
        "diff_prec_raw",
        "diff_prod_avecalib",
        "diff_prod_ave",
        "diff_prod_raw",

        "diff150_prec_avecalib",
        "diff150_prec_ave",
        "diff150_prec_raw",
        "diff150_prod_avecalib",
        "diff150_prod_ave",
        "diff150_prod_raw",


        "ave_prec_avecalib",
        "ave_prec_ave",
        "ave_prec_raw",
        "ave_prod_avecalib",
        "ave_prod_ave",
        "ave_prod_raw",

        "stdev_prec_avecalib",
        "stdev_prec_ave",
        "stdev_prec_raw",
        "stdev_prod_avecalib",
        "stdev_prod_ave",
        "stdev_prod_raw",


    ]


    print("Summerizing for", datasets, "keyword =", keyword)


    counterall = {}

    for dataset in datasets:
        print(dataset)
        counter = {}
        for target in targets:
            counter[target] = 0

        logpath = "../data/{}/log/*".format(dataset)
        paths = glob.glob(logpath)
        for file in paths:
            if not "averaging" in file:
                continue
            with open(file) as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if len(row) < 2:
                        continue
                    if row[1] == "Averaged spectra":
                        counter["spec_averaged"] = int(row[2])
                    if row[1] == "Total raw spectra":
                        counter["spec_raw"] = int(row[2])
                    if row[1] == "Used raw spectra":
                        counter["spec_used"] = int(row[2])
                    if row[1] == "Peaks in averaged MS2 spectra":
                        counter["peak_averaged"] = int(row[2])
                    if row[1] == "Total peaks in raw MS2 spectra":
                        counter["peak_raw"] = int(row[2])
                    if row[1] == "Used peaks in raw MS2 spectra":
                        counter["peak_used"] = int(row[2])
                    """
                    summary	Averaged spectra	2873
                    summary	Peaks in averaged MS2 spectra	32334
                    summary	Peaks in averaged MS1 spectra	2732160
                    summary	Total raw spectra	365333
                    summary	Used raw spectra	339576
                    summary	Total peaks in raw MS2 spectra	8052461
                    summary	Used peaks in raw MS2 spectra	3799297
                    summary	Total peaks in raw MS1 spectra	693279633
                    summary	Used peaks in raw MS1 spectra	292090295
                    """

        counter["ave_peak_used"]= counter["peak_used"]/counter["spec_used"]
        counter["ave_peak_averaged"]= counter["peak_averaged"]/counter["spec_averaged"]
        counter["ave_peak_raw"]= counter["peak_raw"]/counter["spec_raw"]



        #
        # delta mz
        #
        filename = "../data/{}/{}_mzcomparison.csv".format(dataset,dataset)
        paths = glob.glob(logpath)

        prec_avecalib_list = []
        prec_ave_list = []
        prec_raw_list = []

        prod_avecalib_list = []
        prod_ave_list = []
        prod_raw_list = []

        prec150_avecalib_list = []
        prec150_ave_list = []
        prec150_raw_list = []

        prod150_avecalib_list = []
        prod150_ave_list = []
        prod150_raw_list = []

        prec_r_avecalib_list = []
        prec_r_ave_list = []
        prec_r_raw_list = []

        prod_r_avecalib_list = []
        prod_r_ave_list = []
        prod_r_raw_list = []

        with open(filename) as f:
            reader = csv.reader(f, delimiter=",")
            for row in reader:
                if row[0] == "id":
                    continue
                mztemp = float(row[5])
                if row[2] == "p" and row[3] == "averaged":
                    prec_ave_list.append(abs(float(row[13])))
                    prec_avecalib_list.append(abs(float(row[9])))

                    prec_r_ave_list.append(float(row[13]))
                    prec_r_avecalib_list.append(float(row[9]))


                elif row[2] == "p" and row[3] != "averaged":
                    prec_raw_list.append(abs(float(row[12])))

                    prec_r_raw_list.append(float(row[12]))

                elif row[2] == "f" and row[3] == "averaged":
                    prod_ave_list.append(abs(float(row[15])))
                    prod_avecalib_list.append(abs(float(row[11])))

                    prod_r_ave_list.append(float(row[15]))
                    prod_r_avecalib_list.append(float(row[11]))

                elif row[2] == "f" and row[3] != "averaged":
                    prod_raw_list.append(abs(float(row[14])))

                    prod_r_raw_list.append(float(row[14]))
                if mztemp < 150:
                    if row[2] == "p" and row[3] == "averaged":
                        prec150_ave_list.append(abs(float(row[13])))
                        prec150_avecalib_list.append(abs(float(row[9])))

                    elif row[2] == "p" and row[3] != "averaged":
                        prec150_raw_list.append(abs(float(row[12])))

                    elif row[2] == "f" and row[3] == "averaged":
                        prod150_ave_list.append(abs(float(row[15])))
                        prod150_avecalib_list.append(abs(float(row[11])))

                    elif row[2] == "f" and row[3] != "averaged":
                        prod150_raw_list.append(abs(float(row[14])))
             

        counter["diff_prec_avecalib"] = sum(prec_avecalib_list)/len(prec_avecalib_list)
        counter["diff_prec_ave"]= sum(prec_ave_list)/len(prec_ave_list)
        counter["diff_prec_raw"]= sum(prec_raw_list)/len(prec_raw_list)
        counter["diff_prod_avecalib"]= sum(prod_avecalib_list)/len(prod_avecalib_list)
        counter["diff_prod_ave"]= sum(prod_ave_list)/len(prod_ave_list)
        counter["diff_prod_raw"]  = sum(prod_raw_list)/len(prod_raw_list)      

        counter["diff150_prec_avecalib"] = sum(prec150_avecalib_list)/len(prec150_avecalib_list)
        counter["diff150_prec_ave"]= sum(prec150_ave_list)/len(prec150_ave_list)
        counter["diff150_prec_raw"]= sum(prec150_raw_list)/len(prec150_raw_list)
        counter["diff150_prod_avecalib"]= sum(prod150_avecalib_list)/len(prod150_avecalib_list)
        counter["diff150_prod_ave"]= sum(prod150_ave_list)/len(prod150_ave_list)
        counter["diff150_prod_raw"]  = sum(prod150_raw_list)/len(prod150_raw_list)     

        counter["ave_prec_avecalib"] = sum(prec_r_avecalib_list)/len(prec_r_avecalib_list)
        counter["ave_prec_ave"]= sum(prec_r_ave_list)/len(prec_r_ave_list)
        counter["ave_prec_raw"]= sum(prec_r_raw_list)/len(prec_r_raw_list)
        counter["ave_prod_avecalib"]= sum(prod_r_avecalib_list)/len(prod_r_avecalib_list)
        counter["ave_prod_ave"]= sum(prod_r_ave_list)/len(prod_r_ave_list)
        counter["ave_prod_raw"]  = sum(prod_r_raw_list)/len(prod_r_raw_list)     

        counter["stdev_prec_avecalib"] = statistics.stdev(prec_r_avecalib_list)
        counter["stdev_prec_ave"]= statistics.stdev(prec_r_ave_list)
        counter["stdev_prec_raw"]= statistics.stdev(prec_r_raw_list)
        counter["stdev_prod_avecalib"]= statistics.stdev(prod_r_avecalib_list)
        counter["stdev_prod_ave"]= statistics.stdev(prod_r_ave_list)
        counter["stdev_prod_raw"]  = statistics.stdev(prod_r_raw_list) 


        print(counter)

        #
        # filtered paak number
        #
        datapath = ddatoolbox.paths.generate_paths_data(dataset)
        filelist = datapath["step4inputlist"]
        peak_filtered_list = []
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
            # Spectra cleaning up
            #
            production_correlation_alpha = 0.01
            cleanedspectra = mb2.wash_spectra(
                alpha=production_correlation_alpha, intthres=0
            )
            if len(cleanedspectra) >= 3:
                peak_filtered_list.append(len(cleanedspectra))
        counter["ave_peak_filtered"] = sum(peak_filtered_list)/len(peak_filtered_list)
        counter["spec_filtered"] = len(peak_filtered_list)
        counter["peak_filtered"] = sum(peak_filtered_list)

        print(counter)
        counterall[dataset] = counter


    output = []
    outputtemp = ["dataset"]
    for target in targets:
        outputtemp.append(target)

    output.append(outputtemp)
    print(outputtemp)

    for dataset, counter in counterall.items():
        outputtemp = [dataset]
        for target in targets:
            outputtemp.append(counter[target])

        output.append(outputtemp)
        print(outputtemp)
    with open("statisticsprocessing.txt", "w", encoding="UTF-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(output)

    for outputtemp in output:
        ddatoolbox.ddatoolbox.writeLog(logfilename, outputtemp)
        print(outputtemp)

def MsfinderStructureDBforMetFrag():
    """Convert the MS-FINDER structure database into a MetFrag-compatible local DB.

    Reads MS-FINDER's ``MsfinderStructureDB-VS15.esd`` reference table, derives an
    InChI string for every entry from its SMILES (via RDKit), splits its InChIKey
    into its three dash-separated blocks, and writes the result as
    ``MsfinderStructureDB-VS15_forMetFrag.csv`` in the format MetFrag expects for a
    local structure database (``Identifier``, ``InChI``, ``MonoisotopicMass``,
    ``MolecularFormula``, ``InChIKey1``, ``InChIKey2``, ``SMILES``, ``Name``,
    ``InChIKey3``). Entries whose SMILES cannot be parsed by RDKit are skipped.

    Returns:
        None
    """

    from statistics import stdev
    import pubchempy as pcp
    from rdkit import Chem

    def get_pubchemdatafrominchikey(inchikey):
        """Fetch selected PubChem properties for a compound identified by InChIKey.

        Args:
            inchikey (str): InChIKey of the compound to look up.

        Returns:
            dict: PubChem properties (``IUPACName``, ``MolecularFormula``,
                ``MolecularWeight``, ``IsomericSMILES``, ``ExactMass``,
                ``MonoisotopicMass``, ``inchikey``, ``inchi``) for the
                lowest-CID matching compound, or ``{}`` if no match is found.
        """
        pdcompoundall = pcp.get_compounds(inchikey, 'inchikey')
        cids = sorted([int(x.cid) for x in pdcompoundall])
        if len(cids) >= 1:
            properties = ['IUPACName','MolecularFormula', 'MolecularWeight', 'IsomericSMILES', 'ExactMass', 'MonoisotopicMass', 'inchikey', 'inchi']
            pdproperties = pcp.get_properties(properties, cids[0], 'cid')
            pdproperty = pdproperties[0]
        else:
            pdproperty = {}
        return(pdproperty)


    path = "../software/MSFINDER/Resources/MsfinderStructureDB-VS15.esd"
    outputlist = [["Identifier","InChI","MonoisotopicMass","MolecularFormula","InChIKey1","InChIKey2","SMILES","Name","InChIKey3"]]
    with open(path,  encoding='UTF-8') as f:

        reader = csv.reader(f, delimiter='\t')

        for i, row in enumerate(reader):
            if i == 0:
                continue
            Identifier = i
            MonoisotopicMass = row[4]
            MolecularFormula = row[5]
            SMILES = row[6]
            Name = row[0]

            InChIKey = row[1]
            InChIKey1_temp = InChIKey.split("-")

            InChI = ""
            InChIKey1 = InChIKey1_temp[0]
            InChIKey2 = InChIKey1_temp[1]
            InChIKey3 = InChIKey1_temp[2]

            #pubchemdata = get_pubchemdatafrominchikey(InChIKey)
            #if pubchemdata == {}:
            #   continue
            #print(pubchemdata)
            #InChI = pubchemdata['InChI']


            # Build the molecule
            #smiles = "CCO"  # example SMILES notation for a chemical formula
            mol = Chem.MolFromSmiles(SMILES)

            if mol is not None:
                # Get the InChI
                InChI = Chem.MolToInchi(mol)
                print(Identifier, "InChI:", InChI)
            else:
                print(Identifier, "Failed to build a molecule from the SMILES string.")
                continue



            templist = [Identifier, InChI, MonoisotopicMass, MolecularFormula, InChIKey1, InChIKey2, SMILES, Name, InChIKey3]


            outputlist.append(templist)
            #print(templist)
        with open("MsfinderStructureDB-VS15_forMetFrag.csv", mode="w", encoding='UTF-8', newline="") as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerows(outputlist)




#########################################################################
def statistics_annotation(datasets, spec_option="", keyword=""):


    """Aggregate and summarize annotation-quality (curation-level) statistics per dataset.

    For each dataset, reads every step4 MassBank record and, after the standard
    correlation-based spectral cleanup, checks for a chain of increasingly strict
    annotation-quality criteria: whether a molecular formula/compound/class was
    proposed at all, whether it passed the aggregated (Bayesian) ranking, whether it
    qualifies by SIRIUS Zodiac/Confidence score against an instrument-type-specific
    threshold, and whether its Bayesian-FDR-like q-value is below 0.01. Records are
    bucketed into curation levels (Level 3a-d via ``COMMENT: EVIDENCE_INVENTORY``,
    Level 3/class/4/5 via the above criteria), and metrics such as GNPS spectral-
    entropy/cosine similarity and agreement between the entropy-based and
    Bayesian-ranked top hit are also tallied.

    The summary results are printed to the console and saved as a tab-delimited
    text file named ``"curationcheck_s6zipfdist.txt"``.

    Args:
        datasets (list[str]): Dataset identifiers to process.
            If an empty list is provided, all dataset names under the data folder
            are retrieved via ``ddatoolbox.paths.get_alldataset_names()``.
        spec_option (str, optional): Reserved for specifying spectral processing
            options; currently unused. Defaults to ``""``.
        keyword (str, optional): Filter string; only step4 record files whose
            filename contains ``keyword`` are processed. If empty, all files are
            processed. Defaults to ``""``.

    Returns:
        None

    Note:
        This function calls ``pickle.load()`` without importing the ``pickle``
        module in this file, which will raise a ``NameError`` unless ``pickle`` has
        been made available some other way before this function runs.
    """
    targets = [
        "rawspec",  # total
        "usedspec",  # total
        "total",  # total
        "specok",  # pass wash spec
        "onemf",  # pass formula at least one
        "aggmf",  # pass formula aggregated ranking
        "qualmf",  # pass qualifing
        "fdrmf",  # fdr was evaluated
        "passlevel4",  # pass level 4 formula
        "entropy",  # number of entropy < 0.7
        "cosine",  # number of cosine < 0.7
        "onecp",  # pass compound at least one
        "aggcp",  # pass compound aggregated ranking
        "qualcp",  # pass qualifing
        "fdrcp",  # fdr was evaluated
        "enttop",  # number of entropy < 0.7 == top ranking
        "passclass",  # pass compound level 3
        "passlevel3",  # pass compound level 3
        "entorlevel3",  # number of entropy < 0.7 == level 3
        "entlevel3",  # number of entropy < 0.7 == level 3
        "level3hmdb",  # number of level 3 is in HMDB
        "level5",  # number of level 5
        "level4",  # number of level 4
        "levelclass",  # number of level 4
        "level3",  # number of level 4
        "level3d",  # number compound level 3d
        "level3c",  # number compound level 3a
        "level3b",  # number compound level 3a
        "level3a",  # number compound level 3a
        "uniqueinchikey",  # number compound level 3a
    ]

    InChIkeydb = {}
    with open("../script/res/InChIkeydb.pkl", "rb") as f:
        InChIkeydb = pickle.load(f)
    #
    # reading pubchem_checked
    pubchem_checked = ddatoolbox.resources.load_pubchemchecked()


    print("Summerizing for", datasets, "keyword =", keyword)

    counterall = {}
    inchisetall = set()

    for dataset in ddatasets:
        print(dataset)
        counter = {}
        for target in targets:
            counter[target] = 0

        logpath = "../data/{}/log/*".format(dataset)
        paths = glob.glob(logpath)
        for file in paths:
            if not "averaging" in file:
                continue
            with open(file) as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if len(row) < 2:
                        continue
                    if row[1] == "Averaged spectra":
                        averaged = row[2]
                    if row[1] == "Total raw spectra":
                        raw = row[2]
                    if row[1] == "Used raw spectra":
                        used = row[2]

        counter["rawspec"] = int(raw)
        counter["usedspec"] = int(used)

        datapath = ddatoolbox.paths.generate_paths_data(dataset)
        template = ddatoolbox.massbankio.read_RecordTemplate_dataset(dataset)
        INSTRUMENT_TYPE = template["INSTRUMENT_TYPE"]
        filelist = datapath["step4inputlist"]

        inchiset = set()


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
            passlevel3 = False
            passclass = False
            passlevel4 = False

            accession = mb2.get_accession()  # id
            counter["total"] += 1  # total

            #
            # Spectra cleaning up
            #
            production_correlation_alpha = 0.01
            cleanedspectra = mb2.wash_spectra(
                alpha=production_correlation_alpha, intthres=0
            )
            if len(cleanedspectra) < 3:
                continue

            counter["specok"] += 1  # pass wash spec

            anno = mb2.get_proposed(target="formula")

            if 1 in anno["rank"]:
                counter["onemf"] += 1  # pass formula at least one

            methods = set(anno["methodrank"].keys())

            targetset = set(["Bayesian_ranking:biodatabase", "Aggregated_ranking"])
            targetmethodset = methods & targetset
            if len(targetmethodset) == 1:
                targetmethod = list(targetmethodset)[0]
                tempdata = anno["methodrank"][targetmethod]
                if 1 in tempdata:
                    counter["aggmf"] += 1  # pass formula at least one

            if 1 in anno["rank"]:
                if "SIRIUS6_Zodiac" in anno["rank"][1]:
                    score = anno["rank"][1]["SIRIUS6_Zodiac"]["score"]
                    if "TOF" in INSTRUMENT_TYPE:
                        thres = 0.6
                    elif "FT" in INSTRUMENT_TYPE:
                        thres = 0.2
                    else:
                        thres = 1.0
                    if score > thres:
                        counter["qualmf"] += 1  # pass qualifing

            targetset = set(
                ["q_value_like", "q-value_like", "Bayesian_qvalue:biodatabase"]
            )
            targetmethodset = methods & targetset
            if len(targetmethodset) == 1:
                targetmethod = list(targetmethodset)[0]
                tempdata = anno["methodrank"][targetmethod]

                counter["fdrmf"] += 1  # fdr was evaluated

                if 1 in tempdata:
                    if tempdata[1]["score"] < 0.01:
                        counter["passlevel4"] += 1  # pass level 4 formula
                        passlevel4 = True


            anno = mb2.get_proposed(target="compound")

            entinchikey = ""
            cosineinchikey = ""

            methods = set(anno["methodrank"].keys())

            targetset = set(["Entropy:GNPS"])
            targetmethodset = methods & targetset
            if len(targetmethodset) == 1:
                targetmethod = list(targetmethodset)[0]
                tempdata = anno["methodrank"][targetmethod]
                if 1 in tempdata:
                    if tempdata[1]["score"] > 0.7:
                        counter["entropy"] += 1  # number of entropy < 0.7
                        entinchikey = tempdata[1]["tg"]

            targetset = set(["SIMILARITY:GNPS"])
            targetmethodset = methods & targetset
            if len(targetmethodset) == 1:
                targetmethod = list(targetmethodset)[0]
                tempdata = anno["methodrank"][targetmethod]
                if 1 in tempdata:
                    if tempdata[1]["score"] > 0.7:
                        counter["cosine"] += 1  # number of entropy < 0.7
                        cosineinchikey = tempdata[1]["tg"]

            if 1 in anno["rank"]:
                counter["onecp"] += 1  # pass compound at least one

            targetset = set(
                ["Aggregated_ranking:Bayesian", "Bayesian_ranking:biodatabas"]
            )
            targetmethodset = methods & targetset
            if len(targetmethodset) == 1:
                targetmethod = list(targetmethodset)[0]
                tempdata = anno["methodrank"][targetmethod]
                if 1 in tempdata:
                    counter["aggcp"] += 1  # pass compound aggregated ranking

            if 1 in anno["rank"]:
                if "SIRIUS6_Confidence:biodatabase" in anno["rank"][1]:
                    score = anno["rank"][1]["SIRIUS6_Confidence:biodatabase"]["score"]
                    if "TOF" in INSTRUMENT_TYPE:
                        thres = 0.64
                    elif "FT" in INSTRUMENT_TYPE:
                        thres = 0.6
                    else:
                        thres = 1.0
                    if score > thres:
                        counter["qualcp"] += 1  # pass qualifing

            targetset = set(["q-value_like:Bayesian", "Bayesian_qvalue:biodatabase"])
            targetmethodset = methods & targetset
            if len(targetmethodset) == 1:
                targetmethod = list(targetmethodset)[0]
                tempdata = anno["methodrank"][targetmethod]
                if 1 in tempdata:
                    topinchikey = tempdata[1]["tg"]
                    if topinchikey == entinchikey:
                        counter["enttop"] += 1  # number of entropy < 0.7 == top ranking

                    counter["fdrcp"] += 1  # fdr was evaluated

                    if tempdata[1]["score"] < 0.01:
                        counter["passlevel3"] += 1  # pass compound level 3

                        inchiset.add(topinchikey)
                        inchisetall.add(topinchikey)

                        passlevel3 = True
                        if topinchikey == entinchikey:
                            counter[
                                "entlevel3"
                            ] += 1  # number of entropy < 0.7 and level 3
                            counter["entorlevel3"] += 1
                        else:
                            if entinchikey != "":
                                counter["entorlevel3"] += 1  #

                        topinchikeylong = topinchikey + "-UHFFFAOYSA-N"
                        if topinchikey in InChIkeydb:
                            counter["level3hmdb"] += 1  # number of level 3 is in HMDB
            
            anno = mb2.get_proposed(target="class")
            methods = set(anno["methodrank"].keys())
            targetset = set(["q-value_like:Bayesian", "Bayesian_qvalue:biodatabase"])
            targetmethodset = methods & targetset
            if len(targetmethodset) == 1:
                targetmethod = list(targetmethodset)[0]
                tempdata = anno["methodrank"][targetmethod]
                if 1 in tempdata:
                    if tempdata[1]["score"] < 0.01:
                        counter["passclass"] += 1  # pass compound level 3
                        passclass = True


            confidenceetext = mb2.get_anything("COMMENT:", "CONFIDENCE")
            confidenceetext = confidenceetext[0]
            evidencetext = mb2.get_anything("COMMENT:", "EVIDENCE_INVENTORY")
            if len(evidencetext) == 0:
                evidencetext = [""]
                #print(accession, evidencetext)
            evidencetext = evidencetext[0]

            if passlevel3:
                counter["level3"] += 1  # number compound level 3d


            elif passclass:
                counter["levelclass"] += 1  # number of level 4
            elif passlevel4:
                counter["level4"] += 1  # number of level 4
            else:
                counter["level5"] += 1  # number of level 5
            if "[IS]" in evidencetext:
                counter["level3d"] += 1  # number compound level 3d
            if "[IS, CX]" in evidencetext:
                counter["level3c"] += 1  # number compound level 3a
            if "[IS, SL]" in evidencetext:
                counter["level3b"] += 1  # number compound level 3a
            if "[IS, SL, CX]" in evidencetext:
                counter["level3a"] += 1  # number compound level 3a
        counter["uniqueinchikey"] = len(inchiset)
        counterall[dataset] = counter

    totalcounter = {}
    for target in targets:
        totalcounter[target] = 0

    for dataset in datasettemp:
        for target in targets:
            totalcounter[target] += counterall[dataset][target]

    totalcounter["uniqueinchikey"] = len(inchisetall)

    counterall["total"] = totalcounter
    output = []
    outputtemp = ["dataset"]
    for target in targets:
        outputtemp.append(target)

    output.append(outputtemp)
    print(outputtemp)

    for dataset, counter in counterall.items():
        outputtemp = [dataset]
        for target in targets:
            outputtemp.append(counter[target])

        output.append(outputtemp)
        print(outputtemp)
    with open("curationcheck_s6zipfdist.txt", "w", encoding="UTF-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(output)

def exactmassgenerator():
    """Enumerate CHNOPS molecular formulas and write exact masses to per-nominal-mass files.

    Brute-force generates candidate molecular formulas over C (0-84), H, N, O, S,
    and P atom counts, filtering each candidate with the classic Seven Golden Rules
    heuristics (element ratio bounds, ring/double-bond-equivalent range, mass-
    dependent element caps, and element-combination probability checks) to keep the
    search space tractable while still covering the plausible formula space. For
    each surviving formula, computes the exact monoisotopic mass and appends a
    tab-separated line (``formula, RDBE, exact mass, C, H, N, O, P, S`` counts) to
    ``./res/exactmass/<nominal_mass>.txt`` (and to the neighboring +-1 nominal-mass
    file when the fractional mass is close to the boundary), building the lookup
    tables used by molecular-formula search.

    Returns:
        None

    Note:
        History of tuning passes (see comments in the source): 070316 initial
        version; 070510 changed to split output files by nominal mass; 070805
        changed to exhaustively cover low atom counts; 220815, 220120, 220124
        later threshold adjustments.
    """
    ## Script that generates a list of exact (precise) molecular masses
    #070316
    #070510 Changed to split output files by nominal mass
    #070805 Changed to exhaustively cover low atom counts
    #220815
    #
    Cmax =  84;
    passname ="./res/exactmass";
    element = {
    "H" : 1.0078250319,
    "C" : 12,
    "O" : 15.994914622,
    "N" : 14.0030740074,
    "S" : 31.972071,
    "P" : 30.97376
    }
    totalcounter = 0;
    Ccounter = 0;


    for Ccounter in range(Cmax+1):
        print(Ccounter);
        #
        # Rule #4 & 5 . Hydrogen/Carbon element ratio check, heteroatom ratio check
        #
        Hmin = int(Ccounter * 0.2)
        Hmax = int(Ccounter * 3.1)
        Nmax = int(Ccounter * 1.3)
        Omax = int(Ccounter * 1.2) + 1 #Add 220120
        Pmax = int(Ccounter * 0.3)
        Smax = int(Ccounter * 0.8)
        if (Hmax < Ccounter*2+5):
            Hmax= int(Ccounter*2+5)
        if (Omax < 5):
            Omax= 5
        if (Smax < 5):
            Smax= 5
        if (Nmax < 3):
            Nmax= 3
        if (Pmax < 1):
            Pmax= 1

        if (Pmax > 3):
            Pmax= 3
        Smax= 3



        for Hcounter in range(Hmin,Hmax + 1):
            print(Ccounter, Hcounter);
            for Ncounter in range(Nmax + 1):
                print(Ccounter, Hcounter, Ncounter);
                for Ocounter in range(Omax + 1):
                    for Scounter in range(Smax + 1):
                        for Pcounter in range(Pmax + 1):

                            sgr_flag = 0;
                            #
                            # Calculation of IHD (Index of Hydrogen Deficiency, i.e. RDBE)
                            #
                            ihd = Ccounter + 1 - Hcounter / 2 + Ncounter / 2 + Pcounter / 2
                            if (ihd  < -2):
                                continue
                            if (ihd  > 20):
                                continue

                            #
                            # Rule #2. LEWIS and SENIOR check, substituted by the IHD check above
                            #
                            if ((ihd - int(ihd)) > 0):
                                #continue
                                pass

                            exactmass = 0
                            exactmass += Ccounter * element["C"]
                            exactmass += Hcounter * element["H"]
                            exactmass += Ocounter * element["O"]
                            exactmass += Ncounter * element["N"]
                            exactmass += Scounter * element["S"]
                            exactmass += Pcounter * element["P"]
                            mf = ""
                            if Ccounter != 0:
                                mf = mf + "C" + str(Ccounter)
                            if Hcounter != 0:
                                mf = mf + "H" + str(Hcounter)
                            if Ncounter != 0:
                                mf = mf + "N" + str(Ncounter)
                            if Ocounter != 0:
                                mf = mf + "O" + str(Ocounter)
                            if Scounter != 0:
                                mf = mf + "S" + str(Scounter)
                            if Pcounter != 0:
                                mf = mf + "P" + str(Pcounter)
                            if (exactmass == 0):
                                exactmass = 1.0078250319
                            #
                            # Rule #1 . restrictions for element numbers
                            #
                            if (exactmass < 500):
                                if (Ccounter > 35):# modified original is 29
                                    continue
                                if (Hcounter > 72):
                                    continue
                                if (Ncounter > 10):
                                    continue
                                if (Ocounter > 18):
                                    continue
                                if (Scounter > 7):
                                    continue
                                if (Pcounter > 4):
                                    continue
                            elif (exactmass < 1200):
                                if (Ccounter > 80):# modified original is 66
                                    continue
                                if (Hcounter > 152):# modified original is 126
                                    continue
                                if (Ncounter > 25):
                                    continue
                                if (Ocounter > 27):
                                    continue
                                if (Scounter > 8):
                                    continue
                                if (Pcounter > 6):
                                    continue
                            elif (exactmass < 1800):# Added 220124, no particular justification
                                if (Ccounter > 99):
                                    continue
                                if (Hcounter > 189):
                                    continue
                                if (Ncounter > 37):
                                    continue
                                if (Ocounter > 39):
                                    continue
                                if (Scounter > 8):
                                    continue
                                if (Pcounter > 6):
                                    continue
                            else:
                                continue
                            #
                            # Rule #6 . element probability check
                            #
                            if ((Ncounter > 1) and (Ocounter > 1) and (Pcounter > 1) and(Scounter > 1)):
                                if (Ncounter >= 10):
                                    continue
                                if (Ocounter >= 20):
                                    continue
                                if (Scounter >= 3):
                                    continue
                                if (Pcounter >= 4):
                                    continue
                            if ((Ocounter > 1) and (Pcounter > 1) and(Scounter > 1)):
                                if (Ocounter >= 14):
                                    continue
                                if (Scounter >= 3):
                                    continue
                                if (Pcounter >= 3):
                                    continue
                            if ((Ncounter > 1) and (Pcounter > 1) and(Scounter > 1)):
                                if (Ncounter >= 4):
                                    continue
                                if (Scounter >= 3):
                                    continue
                                if (Pcounter >= 3):
                                    continue
                            if ((Ncounter > 3) and (Pcounter > 3) and(Scounter > 3)):
                                if (Ncounter >= 11):
                                    continue
                                if (Scounter >= 22):
                                    continue
                                if (Pcounter >= 6):
                                    continue
                            if ((Ncounter > 6) and (Pcounter > 6) and(Scounter > 6)):
                                if (Ncounter >= 19):
                                    continue
                                if (Scounter >= 14):
                                    continue
                                if (Pcounter >= 8):
                                    continue

                            #mass = int(exactmass);
                            mass = int(exactmass/1.0007 + 0.5)

                            filename = passname+"/"+str(mass)+".txt"
                            f = open(filename, 'a',  encoding='UTF-8')
                            temptext = mf+"\t"+str(ihd)+"\t"+str(exactmass)+"\t"+str(Ccounter)+"\t"+str(Hcounter)+"\t"+str(Ncounter)+"\t"+str(Ocounter)+"\t"+str(Pcounter)+"\t"+str(Scounter)+"\n"
                            f.write(temptext)
                            f.close()
                            if int(exactmass+0.05) != mass:
                                filename = passname+"/"+str(mass+1)+".txt"
                                f = open(filename, 'a',  encoding='UTF-8')
                                f.write(temptext)
                                f.close()
                            if int(exactmass-0.05) != mass:
                                filename = passname+"/"+str(mass- 1)+".txt"
                                f = open(filename, 'a',  encoding='UTF-8')
                                f.write(temptext)
                                f.close()
                            totalcounter = totalcounter + 1


    print("Total =" + str(totalcounter))

        
if __name__ == '__main__':

    description = """
    ShinMassBank various tools\n
    Example 1:\n
    python DDApipeline.py HSA001 HSA002 -a -p -anno MS-FINDER MetFragMSFINDER MetFragPubchemLite -r -c -b 20 -t 5\n
    \n
    Example 2: pw.txt is a single line text file describing password for your sirius account\n
    python DDApipeline.py HSA001 HSA002 -anno SIRIUS5 SIRIUS6 < pw.txt\n

    """

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("datasets", nargs='*', default=[])    # positional argument
    parser.add_argument("-m", "--msp", action='store_true', help='Generation of .msp files')    # required option
    parser.add_argument("-p", "--stat_preprocessing", action='store_true', help='statisticsprocessing.txt')    # required option
    parser.add_argument("-d", "--MsfinderStructureDBforMetFrag", action='store_true', help='Convert MS-FINDER DB to MetFrag')    # required option

    parser.add_argument("-spec", "--setupspecdb", action='store_true', help='Setup MassBank-db and GNDBdatabase for searching')    # required option
    parser.add_argument("-e", "--exactmassgenerator", action='store_true', help='Setup MF search database')    # required option
    parser.add_argument("-k", "--keyword", default="", help='Keyword in accession for curation test purpose')    # required option


    args = parser.parse_args()    # Parse command-line arguments

    datasets = args.datasets
    keyword = args.keyword

    if args.exactmassgenerator == True:
        #
        # for each dataset
        #
        exactmassgenerator()
        


    if args.msp == True:
        #
        # for each dataset
        #
        if len(datasets) == 0:
            datasettemp = ddatoolbox.paths.get_alldataset_names()
        else:
            datasettemp = datasets

        generage_msp(datasets)
    if args.stat_preprocessing == True:
        #
        # for each dataset
        #
        if len(datasets) == 0:
            datasettemp = ddatoolbox.paths.get_alldataset_names()
        else:
            datasettemp = datasets

        statistics_preprocessing(datasets, keyword=keyword)


    if args.MsfinderStructureDBforMetFrag == True:


        MsfinderStructureDBforMetFrag()

    if args.setupspecdb == True:
        #pathstem = "./res/"
        filename = './res/spectra/ALL_GNPS_NO_PROPOGATED.json'
        targetpath = "./res/spectra/GNPSnegative/"
        ddatoolbox.massbankserver.readGNPSjsonfile(filename, targetpath, digit = 2, ionization = "ESI", polarity = "NEGATIVE", MSTYPE = "2")
        targetpath = "./res/spectra/GNPSpositive/"
        ddatoolbox.massbankserver.readGNPSjsonfile(filename, targetpath, digit = 2, ionization = "ESI", polarity = "POSITIVE", MSTYPE = "2")
        path = "./res/spectra/MassBank-data/*"
        targetpath = "./res/spectra/MassBankpositive/"
        ddatoolbox.massbankserver.readMassBankfiles(path, targetpath, digit = 2, ionization = "ESI", polarity = "POSITIVE", MSTYPE = "MS2")
        targetpath = "./res/spectra/MassBanknegative/"
        ddatoolbox.massbankserver.readMassBankfiles(path, targetpath, digit = 2, ionization = "ESI", polarity = "NEGATIVE", MSTYPE = "MS2")

